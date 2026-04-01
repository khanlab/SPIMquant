"""Object-level statistics QC: distributions of detected objects.

Plots count, volume/size distribution, and equivalent-radius distribution
for detected objects (plaques, cells, etc.) from segmentation region
properties.  Objects are filtered to the stain specified by the ``stain``
wildcard.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc
    subject = snakemake.wildcards.subject

    df = pd.read_parquet(snakemake.input.regionprops)

    # Filter to the requested stain (the aggregated parquet has a 'stain' column)
    if "stain" in df.columns:
        df = df[df["stain"] == stain].copy()

    n_objects = len(df)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(
        f"Object Statistics QC\n"
        f"Subject: {subject}  |  Stain: {stain}  |  Method: {desc}  |  "
        f"Total objects: {n_objects:,}",
        fontsize=12,
        fontweight="bold",
    )

    # If no objects were detected, show an empty panel
    if n_objects == 0:
        for ax in axes.flat:
            ax.text(
                0.5,
                0.5,
                "No objects detected",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=14,
                color="gray",
            )
            ax.axis("off")
        plt.tight_layout()
        plt.savefig(snakemake.output.png, dpi=150, bbox_inches="tight")
        plt.close()
        return

    has_nvoxels = "nvoxels" in df.columns

    # --- Panel 1: object size (nvoxels) linear scale ---
    ax = axes[0, 0]
    if has_nvoxels:
        nvoxels = df["nvoxels"].values.astype(float)
        ax.hist(nvoxels, bins=50, color="steelblue", alpha=0.75, edgecolor="white")
        ax.set_xlabel("Volume (voxels)")
        ax.set_ylabel("Count")
        ax.set_title("Object size distribution")
    else:
        ax.text(
            0.5,
            0.5,
            "nvoxels not available",
            ha="center",
            va="center",
            transform=ax.transAxes,
            color="gray",
        )
        ax.axis("off")

    # --- Panel 2: object size log scale ---
    ax = axes[0, 1]
    if has_nvoxels:
        ax.hist(
            np.log10(nvoxels + 1),
            bins=50,
            color="darkorange",
            alpha=0.75,
            edgecolor="white",
        )
        ax.set_xlabel("log\u2081\u2080(volume in voxels + 1)")
        ax.set_ylabel("Count")
        ax.set_title("Object size distribution (log scale)")
    else:
        ax.axis("off")

    # --- Panel 3: equivalent spherical radius ---
    ax = axes[1, 0]
    if has_nvoxels:
        # Radius distribution.
        # Assumes objects are approximately spherical and voxels are isotropic.
        # Formula: r = (3V / 4π)^(1/3), where V is volume in voxels.
        # For anisotropic voxels this is an approximation; physical radius
        # would require multiplying by the voxel size in each dimension.
        radii = (3.0 * nvoxels / (4.0 * np.pi)) ** (1.0 / 3.0)
        ax.hist(radii, bins=50, color="forestgreen", alpha=0.75, edgecolor="white")
        ax.set_xlabel("Equivalent radius (voxels)")
        ax.set_ylabel("Count")
        ax.set_title("Object radius distribution (spherical approximation)")
    else:
        ax.axis("off")

    # --- Panel 4: summary statistics ---
    ax = axes[1, 1]
    ax.axis("off")
    if has_nvoxels:
        summary_text = (
            f"Object count:     {n_objects:>12,}\n"
            f"Mean volume:      {np.mean(nvoxels):>12.1f} vx\n"
            f"Median volume:    {np.median(nvoxels):>12.1f} vx\n"
            f"Max volume:       {np.max(nvoxels):>12.1f} vx\n"
            f"Total volume:     {np.sum(nvoxels):>12.0f} vx\n"
            f"Mean radius:      {np.mean(radii):>12.2f} vx\n"
            f"Median radius:    {np.median(radii):>12.2f} vx"
        )
    else:
        summary_text = f"Object count: {n_objects:,}"
    ax.text(
        0.1,
        0.55,
        summary_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8),
    )
    ax.set_title("Summary statistics")

    plt.tight_layout()
    plt.savefig(snakemake.output.png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved object statistics QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
