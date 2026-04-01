"""Segmentation / vessel overview QC: whole-brain slice montage with mask overlay.

Generates a multi-panel figure with sample slices in three anatomical
orientations (axial, coronal, sagittal), each with the binary segmentation or
vessel mask overlaid on the SPIM background image, plus a max-intensity
projection (MIP) column for each orientation.

Data are loaded via ZarrNii so that the correct physical resolution and
aspect ratio are applied automatically.  The SPIM is loaded with
``downsample_near_isotropic=True`` to obtain near-isotropic voxels.
Voxel spacings are read with ``get_zooms()`` and used to compute the correct
``imshow`` aspect ratio for each panel.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
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
    """Compute imshow aspect ratio for a slice through *step_axis*.

    ZarrNii data is indexed (x, y, z) and ``np.rot90`` is applied before
    display, so the displayed image rows and columns are:

    - step_axis=2 (Z-slice): rot90 of data[:, :, z] → rows=y, cols=x  → dy/dx
    - step_axis=1 (Y-slice): rot90 of data[:, y, :] → rows=z, cols=x  → dz/dx
    - step_axis=0 (X-slice): rot90 of data[x, :, :] → rows=z, cols=y  → dz/dy
    """
    dx, dy, dz = float(zooms[0]), float(zooms[1]), float(zooms[2])
    if step_axis == 2:
        return dy / dx
    if step_axis == 1:
        return dz / dx
    return dz / dy  # step_axis == 0


def main():
    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc
    subject = snakemake.wildcards.subject

    spim_img = ZarrNii.from_ome_zarr(
        snakemake.input.spim,
        level=snakemake.params.level,
        downsample_near_isotropic=True,
        channel_labels=[stain],
        **snakemake.params.zarrnii_kwargs,
    )
    mask_img = ZarrNii.from_ome_zarr(
        snakemake.input.mask,
        level=snakemake.params.mask_level,
        downsample_near_isotropic=True,
        **snakemake.params.zarrnii_kwargs,
    )

    # Voxel dimensions (mm) for physical aspect-ratio correction
    zooms = spim_img.get_zooms()

    spim_data = spim_img.data[0].compute()  # (X, Y, Z)
    mask_data = mask_img.data[0].compute()  # (X, Y, Z), values 0–100

    # Bring mask to the same grid as SPIM if needed
    if mask_data.shape != spim_data.shape:
        factors = [t / s for t, s in zip(spim_data.shape, mask_data.shape)]
        mask_data = zoom(mask_data, factors, order=1)

    spim_norm = _percentile_norm(spim_data)
    # Mask values are 0–100 (field-fraction percent); normalise to 0–1 for display
    mask_norm = np.clip(mask_data / 100.0, 0.0, 1.0)

    n_slices = 5
    orient_labels = ["Axial (Z)", "Coronal (Y)", "Sagittal (X)"]
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
            spim_sl = spim_norm[tuple(idx)]
            mask_sl = mask_norm[tuple(idx)]

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
        mip_mask = np.rot90(np.max(mask_norm, axis=ax_idx))
        ax.imshow(mip_spim, cmap="gray", vmin=0, vmax=1, aspect=aspect)
        mip_masked = np.ma.masked_where(mip_mask < 0.01, mip_mask)
        ax.imshow(mip_masked, cmap="hot", alpha=0.6, vmin=0, vmax=1, aspect=aspect)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("MIP", fontsize=9)

    plt.savefig(snakemake.output.png, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved segmentation overview QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
