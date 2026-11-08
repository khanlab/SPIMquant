"""Instance-level objects from a binary segmentation mask.

Connected components merge touching objects, so in dense tissue a single
"object" can be a network of hundreds of puncta.  Here each mask voxel is
instead assigned to a seed by a seeded watershed:

1. **Smoothed image** ``g1``: the corrected image blurred with a Gaussian of
   ``sigma_um`` (mask-normalised, as in ``dogthresh``, so the zeroed area
   outside the brain does not bleed in).
2. **Seeds**: local maxima of ``g1`` inside the mask, at least
   ``min_dist_um`` apart, whose prominence (h-maxima dynamic) is at least
   ``h``.  Every connected component of the mask that ends up with no seed
   gets one at its maximum, so no foreground voxel is left unassigned.
3. **Watershed** on ``-g1``, restricted to the mask.

Blocks are processed with overlap; an instance is reported only by the block
whose core contains its seed, so each instance is counted exactly once.

``h`` is in corrected-intensity units.  By default it is derived from the
noise level: the robust sd of the high-pass residual ``g1 - G(bg_sigma_um)(g1)``
in tissue outside the mask is computed per block, and
``h = noise_k * (a low percentile of those sds)``.

None of this depends on how the mask was made, so it applies to any
segmentation method's mask.
"""

import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import h_maxima
from skimage.segmentation import watershed

from dogthresh import normalized_blur, sigma_vox

_MAD_TO_SIGMA = 1.4826


def neighbourhood(radius_um, zooms_um):
    """Ellipsoidal footprint of ``radius_um``, at least 3 voxels wide on every
    axis, so a voxel's face neighbours are always part of its neighbourhood.
    """
    radii = [max(radius_um / z, 1.0) for z in zooms_um]
    half = [int(np.floor(r)) for r in radii]
    grids = np.ogrid[tuple(slice(-h, h + 1) for h in half)]
    return sum((g / r) ** 2 for g, r in zip(grids, radii)) <= 1.0


def smooth(img, sigma_um, zooms_um):
    """Mask-normalised Gaussian blur of ``img`` (tissue = ``img > 0``)."""
    img = np.asarray(img, dtype=np.float32)
    weight = (img > 0).astype(np.float32)
    return normalized_blur(img, weight, sigma_vox(sigma_um, zooms_um))


def detect_seeds(response, mask, zooms_um, min_dist_um, h=None):
    """Local maxima of ``response`` inside ``mask``.

    Maxima closer than ``min_dist_um`` are merged.  With ``h``, only regional
    maxima with a dynamic (height above the saddle to a higher maximum) of at
    least ``h`` are kept.  Returns an ``(n, 3)`` int array of voxel indices.
    """
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return np.zeros((0, 3), dtype=int)
    cand = mask.copy()
    if h is not None:
        floor = float(response[mask].min())
        r = np.where(mask, response, floor).astype(np.float64)
        cand &= h_maxima(r, h).astype(bool)
    footprint = neighbourhood(min_dist_um, zooms_um)
    cand &= response == ndi.maximum_filter(response, footprint=footprint)
    # a flat plateau gives several equal maxima: keep one per plateau
    lab, n = ndi.label(cand, structure=np.ones((3, 3, 3)))
    if n == 0:
        return np.zeros((0, 3), dtype=int)
    idx = ndi.maximum_position(response, lab, index=np.arange(1, n + 1))
    return np.asarray(idx, dtype=int).reshape(-1, 3)


def ensure_seed_per_component(seeds, response, mask):
    """Add the ``response`` maximum of every mask component without a seed."""
    cc, n = ndi.label(mask, structure=np.ones((3, 3, 3)))
    if n == 0:
        return seeds.reshape(-1, 3)
    has = np.zeros(n + 1, dtype=bool)
    if len(seeds):
        has[cc[tuple(seeds.T)]] = True
    missing = np.nonzero(~has[1:])[0] + 1
    if len(missing) == 0:
        return seeds.reshape(-1, 3)
    extra = np.asarray(ndi.maximum_position(response, cc, index=missing), dtype=int)
    return np.concatenate([seeds.reshape(-1, 3), extra.reshape(-1, 3)])


def seeded_watershed(basin, seeds, mask):
    """Labels 1..n in seed order; mask voxels not reached by a seed stay 0."""
    markers = np.zeros(mask.shape, dtype=np.int32)
    if len(seeds):
        markers[tuple(seeds.T)] = np.arange(1, len(seeds) + 1, dtype=np.int32)
    return watershed(basin, markers=markers, mask=mask)


def instance_features(labels, ids, intensity, zooms_um):
    """Per-instance features for the labels in ``ids``.

    Returns a dict of arrays: ``nvoxels``, ``centroid`` (n, 3) in voxel
    indices, ``mean_intensity``, ``max_intensity`` and ``elongation`` (square
    root of the ratio of the largest to smallest eigenvalue of the physical
    second-moment matrix, counting each voxel as a box; 1 for a sphere, and
    for a cube of voxels).
    """
    ids = np.asarray(ids, dtype=np.int64)
    n = len(ids)
    out = {
        "nvoxels": np.zeros(n),
        "centroid": np.zeros((n, 3)),
        "mean_intensity": np.zeros(n),
        "max_intensity": np.zeros(n),
        "elongation": np.ones(n),
    }
    if n == 0:
        return out
    # restrict to voxels of the requested instances, sorted by label
    sel = np.isin(labels, ids)
    zz, yy, xx = np.nonzero(sel)
    lab = labels[zz, yy, xx]
    val = np.asarray(intensity, dtype=np.float64)[zz, yy, xx]
    order = np.argsort(lab, kind="stable")
    lab, val = lab[order], val[order]
    pts = np.stack([zz[order], yy[order], xx[order]], axis=1).astype(np.float64)
    starts = np.searchsorted(lab, ids, side="left")
    ends = np.searchsorted(lab, ids, side="right")
    phys = np.asarray(zooms_um, dtype=np.float64)
    voxel_var = np.diag(phys**2 / 12.0)
    for j, (a, b) in enumerate(zip(starts, ends)):
        if b <= a:
            continue
        p = pts[a:b]
        out["nvoxels"][j] = b - a
        out["centroid"][j] = p.mean(axis=0)
        out["mean_intensity"][j] = val[a:b].mean()
        out["max_intensity"][j] = val[a:b].max()
        if b - a >= 2:
            # second moments of the voxels as boxes, not points, so an
            # instance one voxel thick still has a finite extent on that axis
            cov = np.cov((p * phys).T, bias=True) + voxel_var
            ev = np.linalg.eigvalsh(cov)
            out["elongation"][j] = np.sqrt(ev[-1] / ev[0])
    return out


def block_noise_sd(
    img, mask, sigma_um, bg_sigma_um, zooms_um, n, seed, min_samples=10000
):
    """Robust sd of the high-pass residual ``g1 - G(bg)(g1)`` in tissue outside
    ``mask`` (from up to ``n`` sampled voxels), or NaN if there are too few.

    Computed per block because blocks dense with sub-threshold signal have a
    much larger spread than quiet tissue; a low percentile across blocks
    (``h_from_noise``) then reflects the noise rather than that signal.
    """
    img = np.asarray(img, dtype=np.float32)
    tissue = img > 0
    g1 = smooth(img, sigma_um, zooms_um)
    bg = normalized_blur(
        g1, tissue.astype(np.float32), sigma_vox(bg_sigma_um, zooms_um)
    )
    keep = tissue & ~(np.asarray(mask) > 0)
    # stay away from the brain edge, where the normalised blurs see few voxels
    keep &= ndi.binary_erosion(tissue, iterations=2)
    values = (g1 - bg)[keep]
    if values.size < min_samples:
        return float("nan")
    if values.size > n:
        values = np.random.default_rng(seed).choice(values, n, replace=False)
    values = values.astype(np.float64)
    return _MAD_TO_SIGMA * float(np.median(np.abs(values - np.median(values))))


def h_from_noise(block_sds, noise_k, fallback, percentile=25, min_blocks=4):
    """``h = noise_k * (percentile of the per-block noise sds)``, or
    ``fallback`` if there are too few usable blocks.

    Returns ``(h, noise_sd, source)``.
    """
    sds = np.asarray(block_sds, dtype=np.float64)
    sds = sds[np.isfinite(sds) & (sds > 0)]
    if sds.size < min_blocks:
        return (
            float(fallback),
            float("nan"),
            (
                f"fallback ({sds.size} blocks with enough tissue outside the mask, "
                f"need {min_blocks})"
            ),
        )
    sd = float(np.percentile(sds, percentile))
    return float(noise_k * sd), sd, f"noise (p{percentile:g} of {sds.size} blocks)"


def segment_instances(img, mask, zooms_um, sigma_um, min_dist_um, h):
    """Seeds, labels and smoothed image for one (overlapped) block."""
    mask = np.asarray(mask) > 0
    g1 = smooth(img, sigma_um, zooms_um)
    seeds = detect_seeds(g1, mask, zooms_um, min_dist_um, h=h)
    seeds = ensure_seed_per_component(seeds, g1, mask)
    labels = seeded_watershed(-g1, seeds, mask)
    return seeds, labels, g1


def block_instances(
    img,
    mask,
    core_lo,
    core_hi,
    offset,
    zooms_um,
    sigma_um,
    min_dist_um,
    h,
    min_volume_um3=0.0,
):
    """Instances owned by one block.

    ``img`` and ``mask`` are the block including its overlap; voxels
    ``core_lo:core_hi`` (per axis, block coordinates) are its core, and
    ``offset`` is the global index of the block's first voxel.  Only instances
    whose seed lies in the core are returned, with global voxel coordinates.
    """
    mask = np.asarray(mask) > 0
    core = tuple(slice(a, b) for a, b in zip(core_lo, core_hi))
    if not mask[core].any():
        return None
    img = np.asarray(img, dtype=np.float32)
    seeds, labels, g1 = segment_instances(img, mask, zooms_um, sigma_um, min_dist_um, h)
    if len(seeds) == 0:
        return None
    lo, hi = np.asarray(core_lo), np.asarray(core_hi)
    own = np.all((seeds >= lo) & (seeds < hi), axis=1)
    if not own.any():
        return None
    ids = np.nonzero(own)[0] + 1
    feats = instance_features(labels, ids, img, zooms_um)
    vox_volume = float(np.prod(zooms_um))
    keep = feats["nvoxels"] * vox_volume >= min_volume_um3
    if not keep.any():
        return None
    off = np.asarray(offset, dtype=np.float64)
    owned_seeds = seeds[own][keep]
    return {
        "centroid": feats["centroid"][keep] + off,
        "seed": owned_seeds + off,
        "nvoxels": feats["nvoxels"][keep],
        "mean_intensity": feats["mean_intensity"][keep],
        "max_intensity": feats["max_intensity"][keep],
        "peak_intensity": g1[tuple(owned_seeds.T)].astype(np.float64),
        "elongation": feats["elongation"][keep],
    }


def to_physical(voxel_coords, affine):
    """Voxel indices (array axis order) to physical coordinates, as zarrnii's
    ``compute_region_properties`` does."""
    voxel_coords = np.asarray(voxel_coords, dtype=np.float64).reshape(-1, 3)
    hom = np.column_stack([voxel_coords, np.ones(len(voxel_coords))])
    return (hom @ np.asarray(affine, dtype=np.float64).T)[:, :3]
