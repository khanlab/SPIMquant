"""Initial transform for registering an additional scan to its reference scan.

The additional scan (fixed) often covers only part of the sample, at a higher
resolution, and its physical coordinates need not relate to the reference
scan's (moving), so neither the image centers nor the physical coordinates can
be relied on to start the registration. With ``method=match``, both scans are
resampled to a common grid (the reference scan's axes, at its resolution or
coarser), and the translation that best matches the additional scan within the
reference scan is found by normalized cross-correlation over every position (by
FFT), summed over the channels. The scans are assumed to share their axes (e.g.
the same mount), so only the translation is searched; the rigid and affine
registration that follows refines the rest.

The output is a greedy (RAS) affine that maps points in the fixed image to the
moving image, as greedy -ia expects.
"""

import json

import nibabel as nib
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve


def _load(path):
    img = nib.load(path)
    return np.asarray(img.dataobj, dtype=np.float32), img.affine


def _prep(vol):
    """Log intensities: robust to the contrast differences between objectives."""
    vol = np.log1p(np.clip(vol, 0, None))
    return vol.astype(np.float32)


def _grid_axes(affine):
    """Unit axis directions and spacings of a voxel grid."""
    lin = affine[:3, :3]
    spacing = np.linalg.norm(lin, axis=0)
    return lin / spacing, spacing


def resample_to_grid(vol, affine, target_affine, target_shape, smooth=True):
    """Linearly resample `vol` onto a target grid, smoothing first to avoid aliasing."""
    if smooth:
        src_axes, spacing = _grid_axes(affine)
        tgt_axes, target_spacing = _grid_axes(target_affine)
        # target spacing along each source axis (the axes may be permuted)
        proj = np.abs(src_axes.T @ tgt_axes) @ target_spacing
        sigma = np.maximum(proj / spacing - 1, 0) / 2
        if np.any(sigma > 0):
            vol = ndimage.gaussian_filter(vol, sigma)
    vox2vox = np.linalg.inv(affine) @ target_affine
    return ndimage.affine_transform(
        vol,
        vox2vox[:3, :3],
        offset=vox2vox[:3, 3],
        output_shape=tuple(int(n) for n in target_shape),
        order=1,
        mode="constant",
        cval=0.0,
    )


def common_grid(fixed_affine, fixed_shape, moving_affine, moving_shape, min_res):
    """Grids on the moving image's axes, at a common spacing, covering each image.

    Returns (spacing, fixed grid affine, shape, moving grid affine, shape).
    """
    axes, moving_spacing = _grid_axes(moving_affine)
    spacing = np.maximum(moving_spacing, min_res)

    def grid(affine, shape):
        corners = np.array(
            [
                [i, j, k, 1]
                for i in (0, shape[0])
                for j in (0, shape[1])
                for k in (0, shape[2])
            ],
            dtype=float,
        )
        world = (affine @ corners.T)[:3]
        # coordinates along the grid axes, in mm
        coords = axes.T @ world
        lo, hi = coords.min(axis=1), coords.max(axis=1)
        grid_shape = np.maximum(np.ceil((hi - lo) / spacing).astype(int), 1)
        grid_affine = np.eye(4)
        grid_affine[:3, :3] = axes * spacing
        grid_affine[:3, 3] = axes @ lo
        return grid_affine, grid_shape

    return (
        spacing,
        *grid(fixed_affine, fixed_shape),
        *grid(moving_affine, moving_shape),
    )


def ncc_map(image, template, min_overlap=0.5):
    """Normalized cross-correlation of `template` at every offset within `image`.

    Offsets are in image voxels and may put the template partly outside the image
    (which is zero-padded); entry `o + template.shape - 1` of the result is the
    template's corner placed at image voxel `o`. Only voxels where the template
    is nonzero count, so a template with empty margins still matches, and offsets
    that leave less than `min_overlap` of the template inside the image are
    excluded, as a few voxels at the edge can correlate by chance.
    """
    mask = (template > 0).astype(np.float32)
    n = mask.sum()
    if n < 2:
        raise ValueError("the fixed image is empty")
    t = np.where(mask > 0, template - template[mask > 0].mean(), 0).astype(np.float32)
    t_norm = np.sqrt((t**2).sum())

    flip = (slice(None, None, -1),) * 3
    # sums of image and image^2 under the template mask, at each offset
    s1 = fftconvolve(image, mask[flip], mode="full")
    s2 = fftconvolve(image**2, mask[flip], mode="full")
    num = fftconvolve(image, t[flip], mode="full")
    var = np.clip(s2 - s1**2 / n, 0, None)
    denom = np.sqrt(var) * t_norm
    overlap = fftconvolve(np.ones_like(image), mask[flip], mode="full")
    # positions that see too little of the image to compare are excluded
    valid = (denom > 1e-6 * denom.max()) & (overlap >= min_overlap * n - 0.5)
    return np.where(valid, num / np.where(valid, denom, 1), -np.inf)


def match_translation(fixed_vols, fixed_affine, moving_vols, moving_affine, min_res):
    """Translation (mm, RAS) mapping fixed points to moving points, by NCC search."""
    spacing, f_aff, f_shape, m_aff, m_shape = common_grid(
        fixed_affine, fixed_vols[0].shape, moving_affine, moving_vols[0].shape, min_res
    )
    print(
        f"matching on a common grid of {spacing} mm: fixed {f_shape}, moving {m_shape}"
    )

    score = None
    for fixed, moving in zip(fixed_vols, moving_vols):
        f = resample_to_grid(_prep(fixed), fixed_affine, f_aff, f_shape)
        m = resample_to_grid(_prep(moving), moving_affine, m_aff, m_shape)
        s = ncc_map(m, f)
        score = s if score is None else score + s
    score /= len(fixed_vols)

    best = np.unravel_index(np.argmax(score), score.shape)
    offset = np.array(best) - (np.array(f_shape) - 1)
    print(
        f"best mean NCC {score[best]:.3f} with the fixed image at moving voxel {offset}"
    )

    # fixed grid voxel 0 lands on moving grid voxel `offset`
    p_fixed = f_aff[:3, 3]
    p_moving = (m_aff @ np.append(offset, 1.0))[:3]
    return p_moving - p_fixed, float(score[best])


def center_translation(fixed_affine, fixed_shape, moving_affine, moving_shape):
    def center(affine, shape):
        return (affine @ np.append((np.array(shape[:3]) - 1) / 2, 1.0))[:3]

    return center(moving_affine, moving_shape) - center(fixed_affine, fixed_shape)


def main(snakemake):
    method = snakemake.params.method
    fixed = [_load(p) for p in snakemake.input.fixed]
    moving = [_load(p) for p in snakemake.input.moving]
    fixed_affine, moving_affine = fixed[0][1], moving[0][1]

    ncc = None
    if method == "match":
        t, ncc = match_translation(
            [v for v, _ in fixed],
            fixed_affine,
            [v for v, _ in moving],
            moving_affine,
            min_res=snakemake.params.min_res_um / 1000.0,
        )
    elif method == "centers":
        t = center_translation(
            fixed_affine, fixed[0][0].shape, moving_affine, moving[0][0].shape
        )
    elif method == "physical":
        t = np.zeros(3)
    else:
        raise ValueError(f"unknown initialization method: {method}")

    xfm = np.eye(4)
    xfm[:3, 3] = t
    print(f"initial translation (mm, fixed to moving): {t}")
    np.savetxt(snakemake.output.xfm_ras, xfm, fmt="%.10g")
    with open(snakemake.output.info, "w") as f:
        json.dump({"method": method, "translation_mm": t.tolist(), "ncc": ncc}, f)


if __name__ == "__main__":
    main(snakemake)  # noqa: F821
