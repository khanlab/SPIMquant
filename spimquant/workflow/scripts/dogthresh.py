"""Difference-of-Gaussians (DoG) blob segmentation with pool exclusion.

Detects small, punctate features (e.g. pS129 alpha-synuclein inclusions) on a
local background, while ignoring large, bright regions such as stain pooling
on the brain surface or in sulci, which defeat global intensity thresholds.

Method string: ``dog+s{s}r{r}k{k}`` (``p`` is the decimal separator, so
``dog+s2p5r4k6`` is s=2.5):

* ``s`` -- sigma of the inner Gaussian, in micrometres.
* ``r`` -- ratio of outer to inner sigma; the outer Gaussian estimates the
  local background at ``r * s`` micrometres.
* ``k`` -- threshold in robust standard deviations:
  ``threshold = median + k * 1.4826 * MAD`` of the DoG response in tissue.

Steps:

1. **DoG response** ``G(s) - G(r*s)``, computed chunkwise with
   ``map_overlap``.  The corrected image is zero outside the brain mask, so
   both blurs are mask-normalised (``G(I*M) / G(M)``, with ``M = I > 0``);
   otherwise the step at the mask edge would show up as a bright rim.
2. **Threshold** from the median and MAD of DoG values sampled from a random
   subset of tissue-containing chunks.
3. **Pool mask.**  On a coarser pyramid level, connected components of
   ``intensity > pool_factor * tissue median`` whose volume exceeds
   ``pool_min_volume_um3`` are pools (bright puncta are far smaller).  The
   pool mask is dilated, upsampled and, at full resolution, restricted to
   voxels whose inner-blurred intensity is itself above the pool threshold,
   then dilated again by ``pool_dilate_um`` to cover the DoG rim around the
   pool.  DoG-positive voxels inside it are dropped.

Outputs a 0/100 mask (as the other methods do, so field fraction can be
computed by local averaging), a QC figure and the low-resolution pool mask.
"""

import re

import numpy as np
from scipy import ndimage as ndi

# Fixed seed so the threshold is reproducible run to run
_RANDOM_STATE = 42

# Gaussian kernels are truncated at this many sigmas (scipy's default); the
# map_overlap depth is derived from it so blocks see the whole kernel.
_TRUNCATE = 4.0

# Below this, G(M) is treated as "no tissue nearby" in the normalised blur.
_MIN_WEIGHT = 1e-3

# MAD -> standard deviation for normally distributed data
_MAD_TO_SIGMA = 1.4826

_NUM = r"(\d+(?:p\d+)?)"


def parse_dog_method(method):
    """Parse ``dog+s{s}r{r}k{k}`` into ``(sigma_um, ratio, k)`` floats."""
    match = re.fullmatch(rf"dog\+s{_NUM}r{_NUM}k{_NUM}", method)
    if match is None:
        raise ValueError(
            f"Invalid dog seg_method '{method}'. Expected dog+s<s>r<r>k<k>, "
            "using 'p' as the decimal separator, e.g. dog+s3r4k6 or dog+s2p5r4k6"
        )
    sigma_um, ratio, k = (float(g.replace("p", ".")) for g in match.groups())
    if sigma_um <= 0 or ratio <= 1:
        raise ValueError(
            f"dog seg_method '{method}': s must be > 0 and r must be > 1 "
            f"(got s={sigma_um}, r={ratio})"
        )
    return sigma_um, ratio, k


def sigma_vox(sigma_um, zooms_um):
    """Per-axis sigma in voxels for a sigma in micrometres."""
    return tuple(float(sigma_um) / float(z) for z in zooms_um)


def ellipsoid(radius_um, zooms_um):
    """Boolean ellipsoidal structuring element of ``radius_um`` micrometres.

    Returns ``None`` when the radius is smaller than half a voxel on every
    axis (i.e. no dilation).
    """
    radii = [radius_um / z for z in zooms_um]
    half = [int(np.floor(r)) for r in radii]
    if all(h == 0 for h in half):
        return None
    grids = np.ogrid[tuple(slice(-h, h + 1) for h in half)]
    dist = sum((g / max(r, 1e-6)) ** 2 for g, r in zip(grids, radii))
    return dist <= 1.0


def normalized_blur(img, weight, sigma):
    """Gaussian blur of ``img`` using only voxels where ``weight`` is set."""
    num = ndi.gaussian_filter(img * weight, sigma, truncate=_TRUNCATE)
    den = ndi.gaussian_filter(weight, sigma, truncate=_TRUNCATE)
    out = np.zeros_like(num)
    np.divide(num, den, out=out, where=den > _MIN_WEIGHT)
    return out


def dog_response(img, sigma1, sigma2):
    """Mask-normalised DoG of a 3D block.

    Returns ``(dog, g1, tissue)``, where ``g1`` is the inner-blurred image
    and ``tissue`` is ``img > 0``.  ``dog`` and ``g1`` are zero outside tissue.
    """
    img = np.asarray(img, dtype=np.float32)
    tissue = img > 0
    weight = tissue.astype(np.float32)
    g1 = normalized_blur(img, weight, sigma1)
    g2 = normalized_blur(img, weight, sigma2)
    dog = np.where(tissue, g1 - g2, 0).astype(np.float32)
    g1 = np.where(tissue, g1, 0).astype(np.float32)
    return dog, g1, tissue


def robust_threshold(values, k):
    """``median + k * 1.4826 * MAD`` of ``values``; returns (thr, med, sigma)."""
    values = np.asarray(values, dtype=np.float64)
    med = float(np.median(values))
    sigma = _MAD_TO_SIGMA * float(np.median(np.abs(values - med)))
    return med + k * sigma, med, sigma


def pool_threshold(tissue_median, pool_factor, log_intensity=False):
    """Intensity above which a voxel counts as bright for pool detection:
    ``pool_factor`` x the tissue median, or, for log intensities
    (``--correction_method log``), the same ratio as an offset,
    ``median + log(pool_factor)``."""
    if log_intensity:
        return tissue_median + float(np.log(pool_factor))
    return pool_factor * tissue_median


def pool_mask_lowres(lowres, pool_thresh, min_volume_um3, zooms_um, dilate_um):
    """Large bright connected components ("pools") on a coarse grid.

    Parameters
    ----------
    lowres : ndarray (3D)
        Corrected intensity at a coarse pyramid level.
    pool_thresh : float
        Intensity above which a voxel counts as bright.
    min_volume_um3 : float
        Bright components larger than this are pools.
    zooms_um : sequence of float
        Voxel size of ``lowres`` in micrometres.
    dilate_um : float
        Dilation applied to the pool mask, so that its upsampled version
        covers the full extent of the pool at full resolution.
    """
    bright = lowres > pool_thresh
    labels, n_labels = ndi.label(bright)
    if n_labels == 0:
        return np.zeros(lowres.shape, dtype=bool)
    vox_vol = float(np.prod(zooms_um))
    sizes = np.bincount(labels.ravel()) * vox_vol
    sizes[0] = 0
    pool = (sizes > min_volume_um3)[labels]
    struct = ellipsoid(dilate_um, zooms_um)
    if struct is not None and pool.any():
        pool = ndi.binary_dilation(pool, structure=struct)
    return pool


def segment_block(img, pool_up, sigma1, sigma2, threshold, pool_thresh, pool_struct):
    """Full-resolution DoG mask (0/100 uint8) for one overlapped block.

    ``pool_up`` is the upsampled low-resolution pool mask for the same block.
    Within it, the pool is refined to voxels whose inner-blurred intensity
    exceeds ``pool_thresh``, then dilated by ``pool_struct``; DoG-positive
    voxels inside are dropped.
    """
    dog, g1, tissue = dog_response(img, sigma1, sigma2)
    mask = (dog > threshold) & tissue
    pool_up = np.asarray(pool_up, dtype=bool)
    if pool_up.any():
        pool = pool_up & (g1 > pool_thresh)
        if pool_struct is not None and pool.any():
            pool = ndi.binary_dilation(pool, structure=pool_struct)
        mask &= ~pool
    return mask.astype(np.uint8) * 100


def lowres_index(block_start, block_shape, full_shape, lowres_shape):
    """Per-axis lowres voxel indices for one full-resolution block.

    Pyramid levels share an origin and differ by a per-axis scale of
    ``full_shape / lowres_shape``, so full-resolution voxel ``i`` maps to
    lowres voxel ``floor((i + 0.5) * lowres_n / full_n)``.
    """
    idx = []
    for start, n_block, n_full, n_low in zip(
        block_start, block_shape, full_shape, lowres_shape
    ):
        i = np.arange(start, start + n_block)
        j = np.floor((i + 0.5) * n_low / n_full).astype(np.int64)
        idx.append(np.clip(j, 0, n_low - 1))
    return idx


def upsample_nearest_block(lowres, block_shape, block_start, full_shape):
    """Nearest-neighbour lookup of ``lowres`` for one full-resolution block.

    ``lowres`` may be any array supporting basic slicing (e.g. a zarr array);
    only the bounding region under the block is read.
    """
    idx = lowres_index(block_start, block_shape, full_shape, lowres.shape)
    region = np.asarray(lowres[tuple(slice(int(i[0]), int(i[-1]) + 1) for i in idx)])
    return region[np.ix_(*(i - i[0] for i in idx))]


def tissue_fraction_per_block(lowres, chunks, full_shape):
    """Fraction of nonzero lowres voxels under each full-resolution chunk.

    Returns a dict mapping block index tuples to tissue fractions.
    """
    edges = []
    for axis_chunks, n_full, n_low in zip(chunks, full_shape, lowres.shape):
        starts = np.concatenate([[0], np.cumsum(axis_chunks)])
        lo = np.floor(starts * n_low / n_full).astype(int)
        lo[-1] = n_low
        edges.append(lo)
    nonzero = lowres > 0
    fractions = {}
    for bz in range(len(chunks[0])):
        for by in range(len(chunks[1])):
            for bx in range(len(chunks[2])):
                sub = nonzero[
                    edges[0][bz] : max(edges[0][bz + 1], edges[0][bz] + 1),
                    edges[1][by] : max(edges[1][by + 1], edges[1][by] + 1),
                    edges[2][bx] : max(edges[2][bx + 1], edges[2][bx] + 1),
                ]
                fractions[(bz, by, bx)] = float(sub.mean()) if sub.size else 0.0
    return fractions


def select_stat_blocks(fractions, n_blocks, min_fraction=0.5, seed=_RANDOM_STATE):
    """Choose up to ``n_blocks`` tissue-rich blocks at random.

    Blocks with at least ``min_fraction`` tissue are preferred; if there are
    too few, the remaining blocks with the most tissue are added.
    """
    rng = np.random.default_rng(seed)
    rich = [b for b, f in fractions.items() if f >= min_fraction]
    rng.shuffle(rich)
    chosen = rich[:n_blocks]
    if len(chosen) < n_blocks:
        rest = sorted(
            (b for b, f in fractions.items() if 0 < f < min_fraction),
            key=lambda b: -fractions[b],
        )
        chosen += rest[: n_blocks - len(chosen)]
    return chosen


def _sample_tissue_dog(img, sigma1, sigma2, n_samples, seed):
    """DoG values at up to ``n_samples`` random tissue voxels of a block."""
    dog, _, tissue = dog_response(img, sigma1, sigma2)
    values = dog[tissue]
    if values.size > n_samples:
        rng = np.random.default_rng(seed)
        values = rng.choice(values, size=n_samples, replace=False)
    return values


def load_pool_lowres(corrected_path, params):
    """Low-resolution image, tissue median and pool mask.

    Reads the pyramid level ``params.pool_level_offset`` levels coarser than
    ``corrected_path`` (or the coarsest available).  The lowres image is GBs
    for a whole brain (~2 Gvox at level +2), so it is read and processed in
    this process with the threaded scheduler: with the distributed scheduler,
    compute() would gather it on a single worker, whose memory limit is only
    the job memory divided by the number of workers.

    Returns a dict with ``lowres`` (float32), ``zooms_um``, ``level``,
    ``znimg`` (the lowres ZarrNii), ``tissue_median``, ``pool_thresh``,
    ``log_intensity``, ``pool`` (bool) and ``pool_frac``.
    """
    import dask
    from zarrnii import ZarrNii

    with dask.config.set(scheduler="threads"):
        level = int(params.pool_level_offset)
        while True:
            try:
                znimg_lo = ZarrNii.from_file(corrected_path, level=level)
                break
            except Exception:
                if level == 0:
                    raise
                level -= 1
        lowres = np.asarray(znimg_lo.data[0].compute(), dtype=np.float32)
    zooms_um = [z * 1000.0 for z in znimg_lo.get_zooms()]
    tissue_lo = lowres[lowres > 0]
    if tissue_lo.size == 0:
        raise ValueError("Corrected image has no nonzero voxels")
    tissue_median = float(np.median(tissue_lo))
    del tissue_lo
    log_intensity = params.correction_method == "log"
    pool_thresh = pool_threshold(
        tissue_median, float(params.pool_factor), log_intensity
    )
    pool = pool_mask_lowres(
        lowres,
        pool_thresh,
        float(params.pool_min_volume_um3),
        zooms_um,
        float(params.pool_lowres_dilate_um),
    )
    pool_frac = float(pool[lowres > 0].mean())
    print(
        f"pool mask at level +{level} (voxel {np.round(zooms_um, 1)} um): "
        f"tissue median {tissue_median:.4g}, pool threshold {pool_thresh:.4g}, "
        f"{100 * pool_frac:.2f}% of tissue"
    )
    return dict(
        lowres=lowres,
        zooms_um=zooms_um,
        level=level,
        znimg=znimg_lo,
        tissue_median=tissue_median,
        pool_thresh=pool_thresh,
        log_intensity=log_intensity,
        pool=pool,
        pool_frac=pool_frac,
    )


def sample_dog_stats(img, lowres, sigma1, sigma2, n_blocks, n_samples):
    """DoG values sampled from tissue-rich blocks of the dask array ``img``.

    Chooses up to ``n_blocks`` blocks of ``img``'s chunk grid (tissue
    fractions from ``lowres``) and samples ``n_samples`` tissue voxels of the
    DoG in each.  The overlap margin is also sampled; it is real neighbouring
    tissue, so this does not bias the statistics.  Returns
    ``(stat_blocks, values)``; the same chunking gives the same blocks, so
    the DoG tuning report reproduces the k units of a dogthresh run.
    """
    import dask
    import dask.array as da

    fractions = tissue_fraction_per_block(lowres, img.chunks, img.shape)
    stat_blocks = select_stat_blocks(fractions, n_blocks)
    if not stat_blocks:
        raise ValueError("No tissue-containing blocks found for DoG statistics")
    stat_depth = tuple(int(np.ceil(_TRUNCATE * s2)) for s2 in sigma2)
    overlapped = da.overlap.overlap(
        img, depth=stat_depth, boundary=dict.fromkeys(range(3), 0)
    )
    samples = [
        dask.delayed(_sample_tissue_dog)(
            overlapped.blocks[b], sigma1, sigma2, n_samples, _RANDOM_STATE + i
        )
        for i, b in enumerate(stat_blocks)
    ]
    return stat_blocks, np.concatenate(dask.compute(*samples))


def _make_qc_figure(values, med, sigma, threshold, sigma_um, ratio, k, pool_info):
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    lo, hi = np.percentile(values, [0.1, 99.9])
    hi = max(hi, threshold * 1.5)
    ax.hist(
        np.clip(values, lo, hi),
        bins=400,
        range=(lo, hi),
        color="steelblue",
        alpha=0.6,
        label="DoG response (tissue sample)",
    )
    ax.set_yscale("log")
    ax.axvline(med, color="k", linewidth=1, label=f"median = {med:.2f}")
    ax.axvline(
        med + sigma,
        color="orange",
        linestyle="--",
        linewidth=1,
        label=f"median + 1.4826·MAD = {med + sigma:.2f}",
    )
    ax.axvline(
        threshold,
        color="red",
        linewidth=1.8,
        label=f"threshold (k={k:g}) = {threshold:.2f}",
    )
    ax.set_xlabel("DoG response (corrected intensity units)")
    ax.set_ylabel("voxels (log)")
    ax.set_title(
        f"DoG blob threshold: σ1={sigma_um:g} µm, σ2={sigma_um * ratio:g} µm, k={k:g}"
    )
    ax.text(
        0.99,
        0.6,
        pool_info,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        family="monospace",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
    )
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    return fig


def main():
    import tempfile
    from contextlib import ExitStack

    import dask
    import dask.array as da
    import zarr
    from dask_setup import get_dask_client
    from zarrnii import ZarrNii

    params = snakemake.params
    sigma_um, ratio, k = parse_dog_method(snakemake.wildcards.desc)
    tmpdir = snakemake.resources.get("tmpdir", None)

    with ExitStack() as stack:
        stack.enter_context(
            get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads)
        )
        znimg = ZarrNii.from_file(snakemake.input.corrected)
        zooms_um = [z * 1000.0 for z in znimg.get_zooms()]  # mm -> um
        sigma1 = sigma_vox(sigma_um, zooms_um)
        sigma2 = sigma_vox(sigma_um * ratio, zooms_um)
        print(
            f"voxel size (um): {np.round(zooms_um, 3)}; "
            f"sigma1 (vox): {np.round(sigma1, 2)}; sigma2 (vox): {np.round(sigma2, 2)}"
        )

        # ---------------------------------------------------------------- #
        # 1. Low-resolution image: tissue median and pool mask
        # ---------------------------------------------------------------- #
        lo = load_pool_lowres(snakemake.input.corrected, params)
        lowres, pool_lo, pool_thresh = lo["lowres"], lo["pool"], lo["pool_thresh"]
        log_intensity = lo["log_intensity"]
        with dask.config.set(scheduler="threads"):
            znimg_pool = lo["znimg"].copy()
            znimg_pool.data = da.from_array(pool_lo.astype(np.uint8)[np.newaxis] * 100)
            znimg_pool.to_nifti(snakemake.output.pool_nii)

        # ---------------------------------------------------------------- #
        # 2. Full-resolution array, rechunked for bounded per-task memory
        # ---------------------------------------------------------------- #
        chunk = int(params.chunk_size)
        if znimg.data.ndim != 4:
            raise ValueError(
                f"Expected a (c, z, y, x) image, got shape {znimg.data.shape}"
            )
        img = znimg.data[0].rechunk((chunk, chunk, chunk))
        full_shape = img.shape
        pool_struct = ellipsoid(float(params.pool_dilate_um), zooms_um)
        pool_dilate_vox = (
            [0, 0, 0] if pool_struct is None else [s // 2 for s in pool_struct.shape]
        )
        depth = tuple(
            int(np.ceil(_TRUNCATE * s2)) + d for s2, d in zip(sigma2, pool_dilate_vox)
        )
        print(f"chunks {img.chunksize}, overlap depth {depth}")

        # ---------------------------------------------------------------- #
        # 3. Threshold from DoG values sampled in tissue-rich blocks
        # ---------------------------------------------------------------- #
        stat_blocks, values = sample_dog_stats(
            img,
            lowres,
            sigma1,
            sigma2,
            int(params.stat_blocks),
            int(params.stat_samples_per_block),
        )
        del lowres
        threshold, med, sigma = robust_threshold(values, k)
        print(
            f"DoG stats from {len(stat_blocks)} blocks ({values.size:,} voxels): "
            f"median {med:.3f}, robust sigma {sigma:.3f}, threshold {threshold:.3f}"
        )

        pool_info = (
            f"tissue median    {lo['tissue_median']:.{3 if log_intensity else 1}f}\n"
            f"pool threshold   {pool_thresh:.{3 if log_intensity else 1}f} "
            f"({params.pool_factor:g}x median"
            f"{', in log space' if log_intensity else ''})\n"
            f"pool min volume  {float(params.pool_min_volume_um3):.3g} um^3\n"
            f"pool fraction    {100 * lo['pool_frac']:.2f}% of tissue "
            f"(level +{lo['level']})"
        )
        fig = _make_qc_figure(
            values, med, sigma, threshold, sigma_um, ratio, k, pool_info
        )
        fig.savefig(snakemake.output.thresholds_png)

        # ---------------------------------------------------------------- #
        # 4. Segment: DoG > threshold, minus refined + dilated pools
        # ---------------------------------------------------------------- #
        # The lowres pool mask (GBs for a whole brain) is stored as a temporary
        # zarr that each task reads its own small region from, rather than
        # being captured in the task graph and shipped to every worker.
        pool_zarr_path = stack.enter_context(
            tempfile.TemporaryDirectory(suffix="_pool.zarr", dir=tmpdir)
        )
        pool_store = zarr.create_array(
            store=pool_zarr_path,
            shape=pool_lo.shape,
            dtype=bool,
            chunks=tuple(min(n, c) for n, c in zip(pool_lo.shape, (64, 256, 256))),
            overwrite=True,
        )
        pool_store[...] = pool_lo
        lowres_shape = pool_lo.shape

        def _pool_block(block, block_info=None):
            start = [loc[0] for loc in block_info[0]["array-location"]]
            pool_store = zarr.open_array(pool_zarr_path, mode="r")
            assert pool_store.shape == lowres_shape
            return upsample_nearest_block(pool_store, block.shape, start, full_shape)

        pool_up = da.map_blocks(
            _pool_block, img, dtype=bool, meta=np.array((), dtype=bool)
        )

        mask = da.map_overlap(
            segment_block,
            img,
            pool_up,
            depth=depth,
            boundary=0,
            dtype=np.uint8,
            meta=np.array((), dtype=np.uint8),
            sigma1=sigma1,
            sigma2=sigma2,
            threshold=threshold,
            pool_thresh=pool_thresh,
            pool_struct=pool_struct,
        )

        znimg_mask = znimg.copy()
        znimg_mask.data = mask[np.newaxis]
        print("segmenting and saving as ome zarr")
        znimg_mask.to_ome_zarr(
            snakemake.output.mask,
            match_scale_factors_from=snakemake.input.corrected,
            **snakemake.config["zarrnii_out_kwargs"],
        )


if __name__ == "__main__":
    main()
