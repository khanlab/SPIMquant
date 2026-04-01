"""Per-channel intensity histogram QC for SPIM data.

Generates linear and log-scale histograms, a cumulative distribution,
saturation/clip fraction, and summary statistics for a single stain channel.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np

from dask_setup import get_dask_client
from zarrnii import ZarrNii


def main():
    stain = snakemake.wildcards.stain
    level = snakemake.params.level
    hist_bins = snakemake.params.hist_bins
    hist_range = snakemake.params.hist_range

    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        znimg = ZarrNii.from_ome_zarr(
            snakemake.input.spim,
            level=level,
            channel_labels=[stain],
            **snakemake.params.zarrnii_kwargs,
        )
        hist_counts, bin_edges = znimg.compute_histogram(
            bins=hist_bins,
            range=hist_range,
        )

    hist_counts = np.asarray(hist_counts, dtype=float)
    bin_edges = np.asarray(bin_edges, dtype=float)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = bin_edges[1] - bin_edges[0]

    total_voxels = hist_counts.sum()
    max_range = hist_range[1]

    # Determine effective display max (last bin with data, with 5 % headroom)
    nonzero_mask = hist_counts > 0
    if nonzero_mask.any():
        disp_max = float(bin_centers[nonzero_mask][-1]) * 1.05
    else:
        disp_max = max_range
    sat_fraction = (
        float(hist_counts[-1]) / total_voxels * 100 if total_voxels > 0 else 0.0
    )

    # Summary statistics derived from histogram
    if total_voxels > 0:
        mean_val = float(np.sum(bin_centers * hist_counts) / total_voxels)
        cumsum_norm = np.cumsum(hist_counts) / total_voxels
        p50_val = float(
            bin_centers[min(np.searchsorted(cumsum_norm, 0.50), len(bin_centers) - 1)]
        )
        p99_val = float(
            bin_centers[min(np.searchsorted(cumsum_norm, 0.99), len(bin_centers) - 1)]
        )
    else:
        mean_val = p50_val = p99_val = 0.0

    # Percentile-based display bounds for the linear-scale histogram panel
    # X: cap at the 99th percentile value (+ 5 % headroom) to avoid long empty tails
    lin_xlim = p99_val * 1.05 if total_voxels > 0 else max_range
    # Y: cap at the tallest bar in the visible x range (+ 5 % headroom) so the
    # body of the distribution is visible rather than dominated by a background spike
    visible = hist_counts[bin_centers <= lin_xlim]
    lin_ylim = (
        float(visible.max()) * 1.05 if visible.size and visible.max() > 0 else 1.0
    )

    subject = snakemake.wildcards.subject

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(
        f"Intensity Histogram QC\nSubject: {subject}  |  Stain: {stain}",
        fontsize=13,
        fontweight="bold",
    )

    # Panel 1: linear-scale histogram
    ax = axes[0, 0]
    ax.bar(bin_centers, hist_counts, width=bin_width, color="steelblue", alpha=0.75)
    ax.set_xlabel("Intensity")
    ax.set_ylabel("Voxel count")
    ax.set_title("Linear-scale histogram")
    ax.set_xlim(0, lin_xlim)
    ax.set_ylim(0, lin_ylim)

    # Panel 2: log-scale histogram
    ax = axes[0, 1]
    log_counts = np.where(hist_counts > 0, np.log10(hist_counts), np.nan)
    ax.bar(bin_centers, log_counts, width=bin_width, color="darkorange", alpha=0.75)
    ax.set_xlabel("Intensity")
    ax.set_ylabel("log\u2081\u2080(voxel count)")
    ax.set_title("Log-scale histogram")
    ax.set_xlim(0, disp_max)

    # Panel 3: cumulative distribution
    ax = axes[1, 0]
    if total_voxels > 0:
        cumsum_pct = cumsum_norm * 100
        ax.plot(bin_centers, cumsum_pct, color="forestgreen", lw=1.5)
        ax.axvline(
            x=p50_val,
            color="purple",
            linestyle="--",
            alpha=0.7,
            label=f"Median ({p50_val:.1f})",
        )
        ax.axvline(
            x=p99_val,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"99th pctile ({p99_val:.1f})",
        )
        ax.legend(fontsize=8)
    ax.set_xlabel("Intensity")
    ax.set_ylabel("Cumulative voxels (%)")
    ax.set_title("Cumulative distribution")
    ax.set_ylim(0, 105)
    ax.set_xlim(0, disp_max)

    # Panel 4: summary statistics
    ax = axes[1, 1]
    ax.axis("off")
    summary_text = (
        f"Total voxels:      {int(total_voxels):>14,}\n"
        f"Mean intensity:    {mean_val:>14.2f}\n"
        f"Median (50th):     {p50_val:>14.2f}\n"
        f"99th percentile:   {p99_val:>14.2f}\n"
        f"Max range:         {max_range:>14.1f}\n"
        f"Saturation frac.:  {sat_fraction:>13.3f}%"
    )
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
    print(f"Saved intensity histogram QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
