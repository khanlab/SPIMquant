"""Segmentation overview QC: slice montage with field-fraction overlay.

Generates a multi-panel figure with sample slices in three anatomical
orientations (axial, coronal, sagittal), each with the segmentation
field-fraction overlaid on the SPIM background image, plus a
max-intensity projection (MIP) column for each orientation.

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


def _match_shape(source, target_shape):
    """Zoom *source* array to *target_shape* if shapes differ."""
    if source.shape == target_shape:
        return source
    factors = [t / s for t, s in zip(target_shape, source.shape)]
    return zoom(source, factors, order=1)


def main():
    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc
    subject = snakemake.wildcards.subject

    spim_data = nib.load(snakemake.input.spim).get_fdata()
    ff_data = nib.load(snakemake.input.fieldfrac).get_fdata()

    # Bring field-fraction image to the same voxel grid as SPIM
    ff_data = _match_shape(ff_data, spim_data.shape)

    spim_norm = _percentile_norm(spim_data)
    # Field fraction values are 0–100; normalise to 0–1 for display
    ff_norm = np.clip(ff_data / 100.0, 0.0, 1.0)

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
        slice_indices = _sample_slices(spim_data.shape[ax_idx], n_slices)

        for col, sl in enumerate(slice_indices):
            idx = [slice(None)] * 3
            idx[ax_idx] = int(sl)
            spim_sl = np.rot90(spim_norm[tuple(idx)])
            ff_sl = np.rot90(ff_norm[tuple(idx)])

            ax = axes[row, col]
            ax.imshow(spim_sl, cmap="gray", vmin=0, vmax=1, aspect="auto")
            ff_masked = np.ma.masked_where(ff_sl < 0.01, ff_sl)
            ax.imshow(ff_masked, cmap="hot", alpha=0.6, vmin=0, vmax=1, aspect="auto")
            ax.set_xticks([])
            ax.set_yticks([])
            if col == 0:
                ax.set_ylabel(orient_name, fontsize=9)
            ax.set_title(f"sl {sl}", fontsize=7)

        # MIP column (last column)
        ax = axes[row, n_slices]
        mip_spim = np.rot90(np.max(spim_norm, axis=ax_idx))
        mip_ff = np.rot90(np.max(ff_norm, axis=ax_idx))
        ax.imshow(mip_spim, cmap="gray", vmin=0, vmax=1, aspect="auto")
        mip_ff_masked = np.ma.masked_where(mip_ff < 0.01, mip_ff)
        ax.imshow(mip_ff_masked, cmap="hot", alpha=0.6, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("MIP", fontsize=9)

    plt.savefig(snakemake.output.png, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved segmentation overview QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
