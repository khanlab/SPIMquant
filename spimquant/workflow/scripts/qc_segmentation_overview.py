"""Segmentation / vessel overview QC: whole-brain slice montage with mask overlay.

Generates a multi-panel figure with sample slices in three anatomical
orientations (axial, coronal, sagittal), each with the binary segmentation or
vessel mask overlaid on the SPIM background image, plus a max-intensity
projection (MIP) column for each orientation.

Data are loaded via ZarrNii (the SPIM with ``downsample_near_isotropic=True``),
converted to NIfTI and reoriented to RAS with the image affine, so that the
rows show true axial, coronal and sagittal slices whatever the storage order
(ZYX for OME-Zarr and Imaris) and orientation of the scan.  Every panel is
displayed with ``np.rot90`` (superior or anterior up) and the aspect ratio of
the voxel spacings.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from scipy.ndimage import zoom
from zarrnii import ZarrNii


def _percentile_norm(arr, pct_low=1, pct_high=99):
    """Percentile-normalise *arr* to the range [0, 1]."""
    lo = np.percentile(arr, pct_low)
    hi = np.percentile(arr, pct_high)
    if hi > lo:
        return np.clip((arr.astype(float) - lo) / (hi - lo), 0.0, 1.0)
    return np.zeros_like(arr, dtype=float)


def _sample_slices(size, n=5):
    """Return *n* evenly-spaced slice indices in the central 60 % of *size*."""
    start = int(size * 0.2)
    stop = int(size * 0.8)
    return np.linspace(start, stop - 1, n, dtype=int)


def _slice_aspect(zooms, step_axis):
    """imshow aspect ratio for a slice through *step_axis* of a RAS volume.

    The volume is indexed (x, y, z) and every panel is ``np.rot90`` of the
    slice, so the displayed rows and columns are:

    - step_axis=2 (axial): rot90 of data[:, :, z] -> rows=y, cols=x -> dy/dx
    - step_axis=1 (coronal): rot90 of data[:, y, :] -> rows=z, cols=x -> dz/dx
    - step_axis=0 (sagittal): rot90 of data[x, :, :] -> rows=z, cols=y -> dz/dy
    """
    dx, dy, dz = float(zooms[0]), float(zooms[1]), float(zooms[2])
    if step_axis == 2:
        return dy / dx
    if step_axis == 1:
        return dz / dx
    return dz / dy  # step_axis == 0


def _ras_volume(img):
    """(x, y, z) array in RAS orientation and its voxel spacings (mm)."""
    nii = nib.as_closest_canonical(img.to_nifti())
    data = np.asanyarray(nii.dataobj)
    while data.ndim > 3 and data.shape[-1] == 1:
        data = data[..., 0]
    if data.ndim != 3:
        raise ValueError(f"expected a single-channel volume, got shape {data.shape}")
    return data, nii.header.get_zooms()[:3]


def main():
    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc
    subject = snakemake.wildcards.subject

    spim_img = ZarrNii.from_file(
        snakemake.input.spim,
        level=snakemake.params.level,
        downsample_near_isotropic=True,
        channel_labels=[stain],
        **snakemake.params.zarrnii_kwargs,
    )
    mask_img = ZarrNii.from_file(
        snakemake.input.mask,
        level=snakemake.params.mask_level,
        downsample_near_isotropic=True,
    )

    # RAS (x, y, z) arrays, whatever the storage order and scan orientation
    spim_data, zooms = _ras_volume(spim_img)
    mask_data, _ = _ras_volume(mask_img)  # values 0-100

    # Bring mask to the same grid as SPIM if needed
    if mask_data.shape != spim_data.shape:
        factors = [t / s for t, s in zip(spim_data.shape, mask_data.shape)]
        mask_data = zoom(mask_data, factors, order=1)

    spim_norm = _percentile_norm(spim_data)
    # Mask values are 0-100 (percent of each voxel covered at full resolution)
    relative = getattr(snakemake.params, "overlay_scale", "fixed") == "relative"
    if relative:
        # Dense structures such as vessels cover a few percent of nearly every
        # voxel at this level, so a fixed 0-100 scale paints the whole brain.
        # Show where the coverage is high relative to this brain instead.
        tissue = mask_data[mask_data > 0]
        lo, hi = np.percentile(tissue, [75, 99.5]) if tissue.size else (0, 100)
        mask_norm = np.clip((mask_data - lo) / max(hi - lo, 1e-6), 0.0, 1.0)
    else:
        mask_norm = np.clip(mask_data / 100.0, 0.0, 1.0)

    n_slices = 5
    orient_labels = ["Axial", "Coronal", "Sagittal"]
    # axis along which we step through slices: 2 = Z, 1 = Y, 0 = X
    step_axes = [2, 1, 0]

    fig, axes = plt.subplots(3, n_slices + 1, figsize=(18, 9), constrained_layout=True)
    fig.suptitle(
        f"Segmentation Overview QC\n"
        f"Subject: {subject}  |  Stain: {stain}  |  Method: {desc}",
        fontsize=12,
        fontweight="bold",
    )

    for row, (orient_name, ax_idx) in enumerate(zip(orient_labels, step_axes)):
        aspect = _slice_aspect(zooms, ax_idx)
        slice_indices = _sample_slices(spim_data.shape[ax_idx], n_slices)

        for col, sl in enumerate(slice_indices):
            idx = [slice(None)] * 3
            idx[ax_idx] = int(sl)
            spim_sl = np.rot90(spim_norm[tuple(idx)])
            mask_sl = np.rot90(mask_norm[tuple(idx)])

            ax = axes[row, col]
            ax.imshow(spim_sl, cmap="gray", vmin=0, vmax=1, aspect=aspect)
            mask_masked = np.ma.masked_where(mask_sl < 0.01, mask_sl)
            ax.imshow(mask_masked, cmap="hot", alpha=0.6, vmin=0, vmax=1, aspect=aspect)
            ax.set_xticks([])
            ax.set_yticks([])
            if col == 0:
                ax.set_ylabel(orient_name, fontsize=9)
            ax.set_title(f"sl {sl}", fontsize=7)

        # MIP column (last column)
        ax = axes[row, n_slices]
        mip_spim = np.rot90(np.max(spim_norm, axis=ax_idx))
        if relative:
            # a max projection of a relative scale saturates: use the mean
            proj = np.mean(mask_data, axis=ax_idx)
            plo, phi = (
                np.percentile(proj[proj > 0], [75, 99.5])
                if (proj > 0).any()
                else (0, 1)
            )
            mip_mask = np.rot90(np.clip((proj - plo) / max(phi - plo, 1e-6), 0.0, 1.0))
        else:
            mip_mask = np.rot90(np.max(mask_norm, axis=ax_idx))
        ax.imshow(mip_spim, cmap="gray", vmin=0, vmax=1, aspect=aspect)
        mip_masked = np.ma.masked_where(mip_mask < 0.01, mip_mask)
        ax.imshow(mip_masked, cmap="hot", alpha=0.6, vmin=0, vmax=1, aspect=aspect)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Mean projection" if relative else "MIP", fontsize=9)

    plt.savefig(snakemake.output.png, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved segmentation overview QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
