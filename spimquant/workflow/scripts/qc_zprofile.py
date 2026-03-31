"""Z-profile QC: per-slice signal intensity and segmented fraction.

Plots the mean signal intensity and mean field fraction (segmented area)
across Z-slices, revealing depth-dependent artefacts, striping, or uneven
staining/illumination.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np


def _per_slice_stats(data, axis=2):
    """Return mean and std of *data* for each index along *axis*."""
    # Compute reduction axes: all axes except the slice axis
    reduce_axes = tuple(i for i in range(data.ndim) if i != axis)
    means = np.mean(data, axis=reduce_axes)
    stds = np.std(data, axis=reduce_axes)
    return means, stds


def _resample_to_length(arr, target_len):
    """Linearly interpolate *arr* to *target_len* samples."""
    if len(arr) == target_len:
        return arr
    src = np.linspace(0, 1, len(arr))
    dst = np.linspace(0, 1, target_len)
    return np.interp(dst, src, arr)


def main():
    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc
    subject = snakemake.wildcards.subject

    spim_data = nib.load(snakemake.input.spim).get_fdata()
    ff_data = nib.load(snakemake.input.fieldfrac).get_fdata()

    # Compute per-z-slice statistics (NIfTI convention: axis 2 = Z)
    spim_mean_z, spim_std_z = _per_slice_stats(spim_data, axis=2)
    ff_mean_z, _ = _per_slice_stats(ff_data, axis=2)

    n_z = len(spim_mean_z)
    z_idx = np.arange(n_z)

    # Resample field-fraction profile to match SPIM Z length if needed
    ff_mean_z = _resample_to_length(ff_mean_z, n_z)

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle(
        f"Z-Profile QC\nSubject: {subject}  |  Stain: {stain}  |  Method: {desc}",
        fontsize=12,
        fontweight="bold",
    )

    # Panel 1: mean intensity per slice
    ax = axes[0]
    ax.plot(z_idx, spim_mean_z, color="steelblue", lw=1.5, label="Mean intensity")
    ax.fill_between(
        z_idx,
        np.maximum(spim_mean_z - spim_std_z, 0),
        spim_mean_z + spim_std_z,
        alpha=0.25,
        color="steelblue",
        label="\u00b1 1 SD",
    )
    ax.set_ylabel("Mean intensity (a.u.)")
    ax.set_title("Mean signal intensity per Z-slice")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel 2: mean field fraction per slice
    ax = axes[1]
    ax.plot(z_idx, ff_mean_z, color="darkorange", lw=1.5, label="Mean field fraction")
    ax.fill_between(z_idx, 0, ff_mean_z, alpha=0.25, color="darkorange")
    ax.set_xlabel("Z-slice index")
    ax.set_ylabel("Mean field fraction (%)")
    ax.set_title("Mean segmented field fraction per Z-slice")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(snakemake.output.png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved Z-profile QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
