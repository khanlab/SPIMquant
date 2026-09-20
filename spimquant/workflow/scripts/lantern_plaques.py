"""5-fold LANTERN ensemble Abeta plaque segmentation.

Runs the published LANTERN ensemble (``apooladi/lantern-ki3-abeta``) over the
Abeta channel and writes the per-fold agreement map -- the fraction of folds
voting foreground, so ``{0, .2, .4, .6, .8, 1}`` for five folds -- as an
OME-Zarr.  Binarization is deliberately left to a downstream rule so the vote
threshold can be changed without re-running inference.

Deviations from ``vesselfm.py``, and why:

- ``ZarrNii.segment()`` cannot be used.  It dispatches through ``da.map_blocks``
  with disjoint blocks, no block position and a hardcoded ``uint8`` output, so it
  can express neither the stride-64 tile overlap the model requires nor
  brain-mask-based block skipping.  ``da.map_overlap`` is driven directly here.

- The model reads RAW intensities.  It was never trained on N4-corrected or
  otherwise rescaled data, so this path takes ``input.spim`` directly rather than
  the ``desc-corrected*`` store the GMM/Otsu rules consume.

Inference contract (from the model card, must match training exactly):
  tiles of 128^3 at stride 64, per-tile 0.5/99.5 percentile clip then z-score,
  softmax over 2 classes, foreground where class-1 probability >= 0.5 per fold.
  Overlapping tiles are combined by ``merge_mode``: ``max`` preserves the
  original LANTERN behavior, while ``average`` and ``gaussian`` provide more
  conventional sliding-window averaging strategies that can reduce edge-driven
  false positives.
"""

import contextlib
import queue

import dask.array as da
import numpy as np
import torch
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii

from dask_setup import get_dask_client
from lantern_model import ResEncLUNet


def tile_origins(extent, tile, stride):
    """Tile origins covering [0, extent), with the last pulled back to hit the edge."""
    if extent <= tile:
        return [0]
    origins = list(range(0, extent - tile + 1, stride))
    if origins[-1] + tile < extent:
        origins.append(extent - tile)
    return origins


def normalize_tiles(x):
    """Per-tile percentile clip then z-score, matching LANTERN's training transform."""
    out = torch.empty_like(x)
    quantiles = torch.tensor([0.005, 0.995], device=x.device, dtype=x.dtype)
    for i in range(x.shape[0]):
        lo, hi = torch.quantile(x[i].reshape(-1), quantiles)
        clipped = torch.clamp(x[i], lo, hi)
        std = clipped.std()
        out[i] = (clipped - clipped.mean()) / (std if std > 1e-8 else 1.0)
    return out


def gaussian_importance_map(tile, sigma_scale=0.125, eps=1e-8):
    """Separable 3D Gaussian importance map for weighting overlapping tiles."""
    if tile < 1:
        raise ValueError(f"tile must be positive, got {tile}")
    center = (tile - 1) / 2.0
    coords = np.arange(tile, dtype=np.float32) - center
    sigma = max(float(tile) * float(sigma_scale), eps)
    kernel_1d = np.exp(-0.5 * (coords / sigma) ** 2).astype(np.float32)
    kernel_1d /= np.max(kernel_1d)
    weights = (
        kernel_1d[:, None, None]
        * kernel_1d[None, :, None]
        * kernel_1d[None, None, :]
    )
    return np.maximum(weights, np.float32(eps))


def merge_tile_votes(
    votes,
    batch_votes,
    batch,
    tile,
    merge_mode,
    n_folds,
    weight_accum=None,
    importance_map=None,
):
    """Merge one batch of per-tile fold votes into the running output volume.

    ``merge_mode="max"`` preserves LANTERN's published behavior: each tile first
    sums fold votes internally, then overlapping tiles combine by per-voxel max.
    That makes a plaque clipped at one tile edge recoverable from another tile
    that contains it whole, and it preserves the calibration of
    ``plaque_vote_threshold``.

    ``average`` and ``gaussian`` instead combine per-tile *fractions* of folds
    voting foreground across all overlapping tiles. ``average`` weights every
    tile position uniformly, while ``gaussian`` down-weights tile edges and
    emphasizes tile centers where the model had more surrounding context.
    """
    if merge_mode == "max":
        for j, (z, y, x) in enumerate(batch):
            region = votes[z : z + tile, y : y + tile, x : x + tile]
            np.maximum(region, batch_votes[j], out=region)
        return

    if merge_mode not in {"average", "gaussian"}:
        raise ValueError(f"unsupported merge_mode={merge_mode!r}")
    if weight_accum is None:
        raise ValueError("weight_accum is required for average/gaussian merge")

    weights = 1.0 if merge_mode == "average" else importance_map
    if merge_mode == "gaussian" and importance_map is None:
        raise ValueError("importance_map is required for gaussian merge")

    scale = np.float32(1.0 / float(n_folds))
    for j, (z, y, x) in enumerate(batch):
        region = votes[z : z + tile, y : y + tile, x : x + tile]
        region_weights = weight_accum[z : z + tile, y : y + tile, x : x + tile]
        tile_fraction = batch_votes[j].astype(np.float32) * scale
        region += tile_fraction * weights
        region_weights += weights


def predict_volume(vol, nets, device, tile, stride, batch_size, merge_mode):
    """Tiled 5-fold vote map for one 3D block.

    Returns raw integer fold-vote counts for ``merge_mode="max"`` (matching the
    historical behavior) and float32 vote fractions in ``[0, 1]`` for
    ``average`` and ``gaussian``.

    The caller owns `device` exclusively for the duration of this call (see
    DevicePool), so no additional locking is needed here.
    """
    orig_shape = vol.shape
    pad = [max(0, tile - n) for n in orig_shape]
    if any(pad):
        # Blocks at the volume border can be thinner than one tile; pad, then crop back.
        vol = np.pad(vol, [(0, p) for p in pad], mode="edge")

    if merge_mode == "max":
        votes = np.zeros(vol.shape, dtype=np.uint8)
        weight_accum = None
        importance_map = None
    elif merge_mode in {"average", "gaussian"}:
        votes = np.zeros(vol.shape, dtype=np.float32)
        weight_accum = np.zeros(vol.shape, dtype=np.float32)
        importance_map = (
            gaussian_importance_map(tile) if merge_mode == "gaussian" else None
        )
    else:
        raise ValueError(f"unsupported merge_mode={merge_mode!r}")

    n_folds = len(nets)
    coords = [
        (z, y, x)
        for z in tile_origins(vol.shape[0], tile, stride)
        for y in tile_origins(vol.shape[1], tile, stride)
        for x in tile_origins(vol.shape[2], tile, stride)
    ]

    for start in range(0, len(coords), batch_size):
        batch = coords[start : start + batch_size]
        arr = np.stack(
            [vol[z : z + tile, y : y + tile, x : x + tile] for z, y, x in batch]
        )
        tensor = torch.from_numpy(arr).to(device)
        tensor = normalize_tiles(tensor)[:, None]
        batch_votes = torch.zeros(
            (len(batch), tile, tile, tile), dtype=torch.uint8, device=device
        )
        with (
            torch.no_grad(),
            torch.autocast(
                device.type, dtype=torch.bfloat16, enabled=device.type == "cuda"
            ),
        ):
            for net in nets:
                prob = torch.softmax(net(tensor).float(), dim=1)[:, 1]
                batch_votes += (prob >= 0.5).to(torch.uint8)
        batch_votes = batch_votes.cpu().numpy()

        merge_tile_votes(
            votes,
            batch_votes,
            batch,
            tile,
            merge_mode,
            n_folds,
            weight_accum=weight_accum,
            importance_map=importance_map,
        )

    if merge_mode == "max":
        return votes[: orig_shape[0], : orig_shape[1], : orig_shape[2]]
    return (
        votes[: orig_shape[0], : orig_shape[1], : orig_shape[2]]
        / np.maximum(
            weight_accum[: orig_shape[0], : orig_shape[1], : orig_shape[2]], 1e-8
        )
    ).astype(np.float32)


def load_folds(model_paths, device):
    """Instantiate every fold on `device`, once, in the main process.

    ~2 GB of fp32 weights for a 5-fold ensemble, so replicating it per GPU is
    cheap.  DevicePool then guarantees only one dask block uses a given device at
    a time, which bounds peak VRAM without needing a lock inside the forward pass.

    ``weights_only=True`` because these checkpoints are downloaded over the
    network; they carry only tensors and plain dicts, so nothing is lost.
    """
    nets = []
    for path in model_paths:
        ckpt = torch.load(path, map_location="cpu", weights_only=True)
        if ckpt.get("model_type") not in (None, "ResEncLUNet"):
            raise ValueError(
                f"{path} declares model_type={ckpt['model_type']!r}, "
                "but this script only builds ResEncLUNet"
            )
        net = ResEncLUNet(num_classes=2, num_input_channels=1)
        net.load_state_dict(ckpt["state_dict"])
        nets.append(net.eval().to(device))
    return nets


def resolve_devices(requested):
    """The devices to run on, clamped to what is actually visible."""
    if not torch.cuda.is_available():
        print("WARNING: no CUDA device visible; falling back to CPU", flush=True)
        return [torch.device("cpu")]
    available = torch.cuda.device_count()
    if requested > available:
        print(
            f"WARNING: {requested} GPUs requested but only {available} visible; "
            f"using {available}",
            flush=True,
        )
    return [torch.device(f"cuda:{i}") for i in range(max(1, min(requested, available)))]


class DevicePool:
    """Lends out one GPU at a time, with the full ensemble resident on each.

    The ensemble is replicated per device (~2 GB of fp32 weights for 5 folds), so
    each GPU works independently on a whole dask block.  Handing out devices from
    a queue rather than striping blocks by index balances the load for free:
    blocks vary enormously in cost because most of a frame is background, so a
    static assignment would leave one card idle -- the same imbalance LANTERN's
    own two-GPU split had to correct for by splitting on cumulative tile count.
    """

    def __init__(self, devices, model_paths):
        self._free = queue.Queue()
        self.nets = {}
        for device in devices:
            self.nets[device] = load_folds(model_paths, device)
            self._free.put(device)

    @contextlib.contextmanager
    def acquire(self):
        device = self._free.get()
        try:
            yield device, self.nets[device]
        finally:
            self._free.put(device)


def build_brain_block_mask(mask_path, image_shape, chunks):
    """Per-block "does this block touch brain" flags, aligned with the image chunking.

    The brain mask lives at a much coarser pyramid level, so rather than
    interpolating it up to the inference grid it is held in memory and indexed by
    each block's array location -- the same approach LANTERN's whole-brain
    inference uses to skip background tiles.
    """
    mask_znimg = ZarrNii.from_nifti(mask_path, axes_order="ZYX")
    mask_arr = np.asarray(mask_znimg.data).squeeze() > 0
    if mask_arr.ndim != 3:
        raise ValueError(f"expected a 3D brain mask, got shape {mask_arr.shape}")

    scales = [m / i for m, i in zip(mask_arr.shape, image_shape[-3:])]

    def block_has_brain(block_info=None):
        shape = block_info[None]["chunk-shape"]
        location = block_info[None]["array-location"]
        index = []
        for axis, (lo, hi) in enumerate(location[-3:]):
            i0 = int(np.floor(lo * scales[axis]))
            i1 = int(np.ceil(hi * scales[axis]))
            i0 = min(max(i0, 0), mask_arr.shape[axis] - 1)
            i1 = min(max(i1, i0 + 1), mask_arr.shape[axis])
            index.append(slice(i0, i1))
        present = bool(mask_arr[index[0], index[1], index[2]].any())
        return np.broadcast_to(np.bool_(present), shape)

    return da.map_blocks(
        block_has_brain,
        chunks=chunks,
        dtype=bool,
        meta=np.array([], dtype=bool),
    )


def plan_chunks(shape, block, halo):
    """Chunk sizes for `shape` at `block`, with no chunk shorter than `halo`.

    A plain rechunk leaves a remainder chunk that can be far shorter than the
    halo (z=1030 at block 512 leaves a 6-voxel tail), and dask refuses to overlap
    a chunk smaller than the depth.  A short tail is absorbed into its
    predecessor instead.  Every chunk start stays a multiple of `block`, so the
    per-block tile grid stays aligned with the global one.
    """
    chunks = [(shape[0],)]
    for extent in shape[1:]:
        if extent <= block:
            chunks.append((extent,))
            continue
        full, remainder = divmod(extent, block)
        if remainder == 0:
            chunks.append((block,) * full)
        elif remainder < halo:
            chunks.append((block,) * (full - 1) + (block + remainder,))
        else:
            chunks.append((block,) * full + (remainder,))
    return tuple(chunks)


def vote_map(data, brain, predict, tile, stride):
    """Lazily map `predict` over `data`, skipping blocks that carry no brain.

    The halo is `tile - stride`, so every voxel in a block's core is covered by
    at least one tile lying wholly inside the block that produced it -- which is
    what makes the stride-64 overlap correct across block boundaries and not just
    within them.  `brain` is halo-expanded alongside `data`, so the skip test is
    conservative: a block is processed if its core *or* its halo touches brain.

    An axis that fits in a single chunk gets no halo: there is no block boundary
    to bridge, and dask rejects a depth larger than the chunk.

    The overlap-mapped result is always float32 so the weighted merge modes keep
    their fractional precision. The legacy ``max`` path still remains exact
    because small integer vote counts are represented exactly in float32.
    """
    halo = tile - stride

    depth = {}
    for axis, axis_chunks in enumerate(data.chunks):
        if axis == 0 or len(axis_chunks) == 1:
            depth[axis] = 0
            continue
        if min(axis_chunks) < halo:
            raise ValueError(
                f"axis {axis} has a chunk of {min(axis_chunks)} voxels, shorter "
                f"than the {halo}-voxel halo; use plan_chunks() to lay out the "
                f"chunks before calling vote_map()"
            )
        depth[axis] = halo

    def segment_block(image_block, brain_block):
        if not brain_block.any():
            return np.zeros(image_block.shape, dtype=np.float32)
        out = np.empty(image_block.shape, dtype=np.float32)
        for c in range(image_block.shape[0]):
            out[c] = predict(image_block[c].astype(np.float32))
        return out

    return da.map_overlap(
        segment_block,
        data,
        brain,
        depth=depth,
        boundary="none",
        trim=True,
        allow_rechunk=False,
        dtype=np.float32,
        meta=np.array([], dtype=np.float32),
    )


def main():
    tile = int(snakemake.params.tile)
    stride = int(snakemake.params.stride)
    merge_mode = str(snakemake.params.merge_mode)
    batch_size = int(snakemake.params.batch_size)
    chunk = int(snakemake.params.chunk)
    n_gpus = int(snakemake.params.n_gpus)

    devices = resolve_devices(n_gpus)

    with get_dask_client("threads", snakemake.threads):
        znimg = ZarrNii.from_file(
            snakemake.input.spim,
            level=int(snakemake.wildcards.level),
            channel_labels=[snakemake.wildcards.stain],
            **snakemake.params.zarrnii_kwargs,
        )

        if znimg.data.ndim != 4 or znimg.data.shape[0] != 1:
            raise ValueError(
                f"expected a single-channel (c,z,y,x) image, got shape {znimg.data.shape}"
            )
        data = znimg.data.rechunk(plan_chunks(znimg.data.shape, chunk, tile - stride))

        brain = build_brain_block_mask(snakemake.input.mask, data.shape, data.chunks)

        pool = DevicePool(devices, snakemake.input.models)
        n_folds = len(pool.nets[devices[0]])
        print(
            f"loaded {n_folds} folds on each of {[str(d) for d in devices]}; "
            f"grid {data.shape} chunk {chunk} halo {tile - stride} "
            f"tile {tile} stride {stride} batch {batch_size} merge {merge_mode}",
            flush=True,
        )

        def predict(vol):
            with pool.acquire() as (device, nets):
                return predict_volume(
                    vol, nets, device, tile, stride, batch_size, merge_mode
                )

        votes = vote_map(data, brain, predict, tile, stride)

        # The fraction of folds voting foreground: float32 in [0, 1]. Under the
        # historical `max` merge it remains in {0, .2, .4, .6, .8, 1} for a 5-fold
        # ensemble, matching LANTERN's own whole-brain probmask; under `average`
        # or `gaussian` it is a weighted mean of those per-tile fractions.
        #
        # Deliberately a probability rather than the raw count or a 0-100 rescaling:
        # it is directly readable as model confidence in a viewer, and independent of
        # how many folds the ensemble happens to have. The 0/100 mask convention does
        # not apply -- that exists so fieldfrac can mean-pool a *mask* into a
        # percentage, and nothing but binarize_lantern_plaques reads this store.
        #
        # float32 costs ~4 bytes/voxel, so a level-1 whole brain is ~100 GB raw; it
        # compresses hard, being overwhelmingly zero.
        znimg_prob = znimg.copy()
        if merge_mode == "max":
            znimg_prob.data = votes.astype(np.float32) / float(n_folds)
        else:
            znimg_prob.data = votes.astype(np.float32)

        with ProgressBar():
            znimg_prob.to_ome_zarr(
                snakemake.output.probseg,
                max_layer=5,
                match_scale_factors_from=snakemake.input.spim,
                **snakemake.config["zarrnii_out_kwargs"],
            )


if __name__ == "__main__":
    main()
