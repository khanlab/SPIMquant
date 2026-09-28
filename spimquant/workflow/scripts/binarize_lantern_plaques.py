"""Threshold the LANTERN vote map, filter components, upsample to the segmentation level.

The ensemble runs near ``plaque_iso_res`` um isotropic because it is
scale-sensitive and was trained at that scale; the segmentation level output is
therefore the upsampled prediction, NOT native inference at that level.

Only the plaque channel of the probseg reaches the mask.  Two optional filters
are decided on the near-iso grid, on connected components (26-connectivity)
of the thresholded plaque channel, before upsampling:

- ``min_size``: components smaller than this many near-iso voxels are dropped
  (the model cards recommend it: +0.06 object F1 on held-out data).
- ``label_filter`` (3-class ensembles only): components touching a class-2
  ("annotated false positive") component of at least ``label_filter_size``
  voxels are dropped.  Class 2 is where at least ``fp_min_votes`` folds put
  their argmax on it (the model cards say any fold).

Both are computed in one pass over the probseg: each 512^3 block is read once
with a one-voxel halo and labelled with scipy, and only the sparse labelled
voxels are kept (plaque is well under 1 percent of the grid).  Labels are then
merged across blocks by matching each block's halo voxels to the core voxels
of the block that owns them, which is exact under 26-connectivity because
every core voxel's whole neighbourhood lies inside its halo'd block.  The
filtered mask is rebuilt from the surviving voxels, so the probseg is never
read a second time and the label image is never materialised.  The near-iso
probseg keeps everything; only the level-0 mask is filtered.

Upsampling is exact replication (``da.repeat``) rather than interpolation, and
the per-axis factors are derived from the shape ratio so this stays correct when
the pyramid downsamples only x and y.

Output is valued 0/100 rather than 0/1, matching ``gmmthresh`` / ``multiotsu`` /
``threshold`` in segmentation.smk, so that mean-pool downsampling in ``fieldfrac``
yields a percentage directly.
"""

import itertools
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import dask
import dask.array as da
import numpy as np
from dask.diagnostics import ProgressBar
from scipy.ndimage import binary_dilation
from scipy.ndimage import label as ndi_label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from zarrnii import ZarrNii

STRUCTURE = np.ones((3, 3, 3), dtype=bool)  # 26-connectivity, as regionprops uses
OFFSETS = [o for o in itertools.product((-1, 0, 1), repeat=3) if o != (0, 0, 0)]


def upsample_nearest(arr, target_shape):
    """Replicate voxels to `target_shape`, then trim/pad to land on it exactly."""
    if len(target_shape) != arr.ndim:
        raise ValueError(f"rank mismatch: {arr.shape} vs target {target_shape}")

    out = arr
    for axis, (have, want) in enumerate(zip(arr.shape, target_shape)):
        factor = max(1, int(round(want / have)))
        if factor > 1:
            out = da.repeat(out, factor, axis=axis)

    out = out[tuple(slice(0, min(h, w)) for h, w in zip(out.shape, target_shape))]
    pad = [(0, w - h) for h, w in zip(out.shape, target_shape)]
    if any(after for _, after in pad):
        out = da.pad(out, pad, mode="edge")
    return out


def block_plan(shape, block):
    """(index, core, outer) for every block: core on the grid, outer = core
    grown by one voxel on each side (clipped), so that every 26-neighbour of a
    core voxel is inside outer."""
    plan = []
    counts = [max(1, -(-n // block)) for n in shape]
    for idx in itertools.product(*(range(c) for c in counts)):
        core = tuple(
            slice(i * block, min((i + 1) * block, n)) for i, n in zip(idx, shape)
        )
        outer = tuple(
            slice(max(0, c.start - 1), min(n, c.stop + 1)) for c, n in zip(core, shape)
        )
        plan.append((idx, core, outer))
    return plan


def _sparse(labels, outer, core):
    """(global voxel keys, local labels) of the labelled voxels in the core of
    `labels` (which covers `outer`), and the same for the halo ring."""
    nz = np.nonzero(labels)
    if nz[0].size == 0:
        e = np.zeros(0, dtype=np.int64)
        return (e, e.astype(np.int32)), (e, e.astype(np.int32))
    lab = labels[nz].astype(np.int32)
    g = [n.astype(np.int64) + o.start for n, o in zip(nz, outer)]
    in_core = np.ones(lab.shape, dtype=bool)
    for gi, c in zip(g, core):
        in_core &= (gi >= c.start) & (gi < c.stop)
    key = (g[0] * _KEY_SHAPE[1] + g[1]) * _KEY_SHAPE[2] + g[2]
    return (key[in_core], lab[in_core]), (key[~in_core], lab[~in_core])


_KEY_SHAPE = None  # (Z, Y, X) of the near-iso grid, for voxel keys


def _contacts(plaque_lab, fp_lab):
    """(plaque local label, fp local label) pairs that are 26-adjacent.

    Sparse: only plaque voxels inside the one-voxel dilation of the class-2
    mask can touch it, so the 26 neighbour lookups are gathered for those
    voxels alone rather than compared over the whole block.
    """
    fp_nz = fp_lab > 0
    if not fp_nz.any():
        return np.zeros((0, 2), dtype=np.int32)
    cand = np.nonzero(binary_dilation(fp_nz, structure=STRUCTURE) & (plaque_lab > 0))
    if cand[0].size == 0:
        return np.zeros((0, 2), dtype=np.int32)
    pl = plaque_lab[cand]
    sh = fp_lab.shape
    pairs = []
    for dz, dy, dx in OFFSETS:
        z, y, x = cand[0] + dz, cand[1] + dy, cand[2] + dx
        ok = (z >= 0) & (z < sh[0]) & (y >= 0) & (y < sh[1]) & (x >= 0) & (x < sh[2])
        fl = fp_lab[z[ok], y[ok], x[ok]]
        m = fl > 0
        if m.any():
            pairs.append(np.stack([pl[ok][m], fl[m]], axis=1))
    if not pairs:
        return np.zeros((0, 2), dtype=np.int32)
    return np.unique(np.concatenate(pairs), axis=0).astype(np.int32)


def analyse_block(votes, item, plaque_thr, fp_thr, want_fp):
    """Read one halo-grown block once, label it, and return only sparse facts."""
    idx, core, outer = item
    v = np.asarray(votes[(slice(None),) + outer])
    plaque = v[0] >= plaque_thr
    out = {"idx": idx}
    if not plaque.any():
        out["plaque"] = None
        return out
    lab, n = ndi_label(plaque, structure=STRUCTURE)
    out["plaque"] = _sparse(lab, outer, core)
    out["n_plaque"] = int(n)
    if want_fp:
        fp = v[1] >= fp_thr
        flab, nf = ndi_label(fp, structure=STRUCTURE)
        out["fp"] = _sparse(flab, outer, core)
        out["n_fp"] = int(nf)
        out["contacts"] = _contacts(lab, flab)
    return out


def merge_labels(blocks, which):
    """Global component id per (block, local label), by joining each block's
    halo voxels to the core voxels of the block that owns them.

    Returns (offsets per block index, component id per global label, and the
    core voxel keys with their component ids), all exact under 26-connectivity
    because every core voxel's full neighbourhood was inside its halo'd block.
    """
    n_key = f"n_{which}"
    offsets = {}
    total = 0
    for b in blocks:
        if b.get(which) is None:
            continue
        offsets[b["idx"]] = total
        total += b[n_key]
    core_keys, core_gid, halo_keys, halo_gid = [], [], [], []
    for b in blocks:
        if b.get(which) is None:
            continue
        (ck, cl), (hk, hl) = b[which]
        off = offsets[b["idx"]]
        core_keys.append(ck)
        core_gid.append(cl.astype(np.int64) + off)
        halo_keys.append(hk)
        halo_gid.append(hl.astype(np.int64) + off)
    core_keys = np.concatenate(core_keys) if core_keys else np.zeros(0, np.int64)
    core_gid = np.concatenate(core_gid) if core_gid else np.zeros(0, np.int64)
    halo_keys = np.concatenate(halo_keys) if halo_keys else np.zeros(0, np.int64)
    halo_gid = np.concatenate(halo_gid) if halo_gid else np.zeros(0, np.int64)
    order = np.argsort(core_keys, kind="stable")
    core_keys, core_gid = core_keys[order], core_gid[order]
    # each halo voxel is a core voxel of exactly one neighbouring block
    pos = np.searchsorted(core_keys, halo_keys)
    pos = np.clip(pos, 0, max(0, core_keys.size - 1))
    ok = core_keys.size > 0
    hit = ok & (core_keys[pos] == halo_keys) if halo_keys.size else np.zeros(0, bool)
    a, b_ = halo_gid[hit], core_gid[pos[hit]]
    n = total + 1
    graph = coo_matrix((np.ones(a.size, dtype=np.int8), (a, b_)), shape=(n, n))
    _, comp = connected_components(graph, directed=False)
    comp[0] = 0
    return offsets, comp, core_keys, comp[core_gid]


def filtered_plaque_mask(
    votes,
    plaque_thr,
    fp_thr,
    min_size,
    label_filter,
    label_filter_size,
    block,
    n_workers,
):
    """Plaque mask on the near-iso grid with the small and the class-2-touching
    components removed, as a lazy dask array built from the sparse survivors."""
    global _KEY_SHAPE
    shape = votes.shape[1:]
    _KEY_SHAPE = tuple(int(n) for n in shape)
    plan = block_plan(shape, block)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        blocks = list(
            pool.map(
                lambda item: analyse_block(
                    votes, item, plaque_thr, fp_thr, label_filter
                ),
                plan,
            )
        )
    n_read = sum(1 for b in blocks if b.get("plaque") is not None)
    print(
        f"  read + labelled {len(plan)} blocks ({n_read} with plaque) in {time.time() - t0:.0f} s",
        flush=True,
    )

    offsets, comp, core_keys, core_comp = merge_labels(blocks, "plaque")
    sizes = np.bincount(core_comp, minlength=comp.max() + 1)
    n_comp = int((sizes[1:] > 0).sum())
    drop = np.zeros(sizes.size, dtype=bool)
    print(
        f"  plaque components: {n_comp}, voxels {core_keys.size} ({time.time() - t0:.0f} s)",
        flush=True,
    )
    if min_size > 0:
        small = (sizes > 0) & (sizes < min_size)
        small[0] = False
        drop |= small
        print(
            f"  min_size={min_size}: dropping {int(small.sum())} components ({time.time() - t0:.0f} s)",
            flush=True,
        )
    if label_filter:
        f_off, f_comp, _, f_core_comp = merge_labels(blocks, "fp")
        f_sizes = np.bincount(f_core_comp, minlength=f_comp.max() + 1)
        big = f_sizes >= label_filter_size
        big[0] = False
        print(
            f"  class-2 components: {int((f_sizes[1:] > 0).sum())}, of which {int(big.sum())} "
            f"have >= {label_filter_size} voxels ({time.time() - t0:.0f} s)",
            flush=True,
        )
        touching = np.zeros(sizes.size, dtype=bool)
        for b in blocks:
            c = b.get("contacts")
            if c is None or c.shape[0] == 0 or b["idx"] not in f_off:
                continue
            pc = comp[c[:, 0].astype(np.int64) + offsets[b["idx"]]]
            fc = f_comp[c[:, 1].astype(np.int64) + f_off[b["idx"]]]
            touching[pc[big[fc]]] = True
        drop |= touching
        print(
            f"  label_filter: dropping {int(touching.sum())} plaque components touching them ({time.time() - t0:.0f} s)",
            flush=True,
        )

    keep = ~drop[core_comp]
    kept_keys = core_keys[keep]
    print(
        f"  keeping {kept_keys.size} of {core_keys.size} plaque voxels ({time.time() - t0:.0f} s)",
        flush=True,
    )

    # Dense blocks on demand from the sorted survivor keys: no second read.
    Z, Y, X = _KEY_SHAPE
    counts = [max(1, -(-n // block)) for n in shape]
    chunks = tuple(
        tuple(min(block, n - i * block) for i in range(c))
        for n, c in zip(shape, counts)
    )

    def fill(block_info=None):
        loc = block_info[None]["array-location"]
        bshape = block_info[None]["chunk-shape"]
        (z0, z1), (y0, y1), (x0, x1) = loc
        lo = (z0 * Y + y0) * X + x0
        hi = ((z1 - 1) * Y + (y1 - 1)) * X + (x1 - 1)
        sel = kept_keys[
            np.searchsorted(kept_keys, lo) : np.searchsorted(
                kept_keys, hi, side="right"
            )
        ]
        out = np.zeros(bshape, dtype=bool)
        if sel.size:
            z = sel // (Y * X)
            r = sel - z * (Y * X)
            y = r // X
            x = r - y * X
            m = (y >= y0) & (y < y1) & (x >= x0) & (x < x1)
            out[z[m] - z0, y[m] - y0, x[m] - x0] = True
        return out

    return da.map_blocks(fill, chunks=chunks, dtype=bool, meta=np.array([], dtype=bool))


def main():
    n_folds = int(snakemake.params.n_folds)
    vote_threshold = int(snakemake.params.vote_threshold)
    min_size = int(snakemake.params.min_size)
    label_filter = bool(snakemake.params.label_filter)
    label_filter_size = int(snakemake.params.label_filter_size)
    fp_min_votes = int(snakemake.params.fp_min_votes)

    if not 1 <= vote_threshold <= n_folds:
        raise ValueError(
            f"plaque_vote_threshold={vote_threshold} is out of range for an "
            f"{n_folds}-fold ensemble; it must be between 1 and {n_folds}"
        )

    dask.config.set(scheduler="threads", num_workers=snakemake.threads)

    probseg = ZarrNii.from_file(snakemake.input.probseg, level=0)
    ref = ZarrNii.from_file(
        snakemake.input.spim,
        level=int(snakemake.params.target_level),
        channel_labels=[snakemake.wildcards.stain],
        **snakemake.params.zarrnii_kwargs,
    )

    votes = probseg.data
    while votes.ndim > 4 and votes.shape[0] == 1:  # (t, c, z, y, x) -> (c, z, y, x)
        votes = votes[0]
    if votes.ndim != 4:
        raise ValueError(
            f"expected a (c, z, y, x) probseg, got shape {probseg.data.shape}"
        )
    n_channels = votes.shape[0]

    # probseg holds votes/n_folds. Compare against the MIDPOINT between adjacent
    # attainable values rather than the exact fraction: 2.5/5 = 0.5 sits safely
    # between 0.4 and 0.6, so no float representation error can flip a voxel at the
    # boundary the way `>= 3/5` could.
    plaque = votes[0] >= (vote_threshold - 0.5) / n_folds
    fp = votes[1] >= (fp_min_votes - 0.5) / n_folds if n_channels >= 2 else None

    if label_filter and fp is None:
        print(
            "--label_filter requested but the probseg has no class-2 channel "
            "(binary ensemble); skipping that filter",
            flush=True,
        )
        label_filter = False

    print(
        f"votes>={vote_threshold} of {n_folds}, {n_channels} channel(s) | "
        f"min_size={min_size} label_filter={label_filter} "
        f"(class-2 size >= {label_filter_size}, votes >= {fp_min_votes}) | "
        f"{tuple(votes.shape[1:])} -> {ref.data.shape}",
        flush=True,
    )

    if min_size > 0 or label_filter:
        plaque = filtered_plaque_mask(
            votes,
            (vote_threshold - 0.5) / n_folds,
            (fp_min_votes - 0.5) / n_folds,
            min_size,
            label_filter,
            label_filter_size,
            block=512,
            n_workers=snakemake.threads,
        )

    upsampled = upsample_nearest(plaque[np.newaxis], ref.data.shape)

    znimg_mask = ref.copy()
    znimg_mask.data = (upsampled * 100).astype(np.uint8)

    with ProgressBar():
        znimg_mask.to_ome_zarr(
            snakemake.output.mask,
            max_layer=5,
            match_scale_factors_from=snakemake.input.spim,
            **snakemake.config["zarrnii_out_kwargs"],
        )


if __name__ == "__main__":
    main()
