"""Compute Multi-Otsu thresholds from a batch-aggregated histogram.

Reads the aggregated histogram ``.npz`` file (containing ``counts`` and
``bin_edges``), optionally filters to a percentile-based intensity range to
exclude background and saturation artefacts, and then applies multi-level Otsu
thresholding.

Outputs:
- A PNG visualisation of the histogram with Otsu threshold positions marked.
- A TSV file listing every threshold index and the corresponding value.

This is a Snakemake script; the ``snakemake`` object is automatically provided
when executed as part of a Snakemake workflow.
"""

import matplotlib
import numpy as np

matplotlib.use("agg")

from zarrnii.analysis import compute_otsu_thresholds


def main():
    data = np.load(snakemake.input.histogram)
    hist_counts = data["counts"].astype(np.float64)
    bin_edges = data["bin_edges"]

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    pct_lo, pct_hi = snakemake.params.hist_percentile_range

    # Derive a percentile-based intensity range from the aggregated histogram
    total = hist_counts.sum()
    if total > 0:
        cumsum_norm = np.cumsum(hist_counts) / total
        idx_lo = int(np.searchsorted(cumsum_norm, pct_lo / 100.0))
        idx_hi = int(np.searchsorted(cumsum_norm, pct_hi / 100.0))
        idx_lo = max(0, min(idx_lo, len(bin_centers) - 1))
        idx_hi = max(idx_lo + 1, min(idx_hi, len(bin_centers) - 1))
    else:
        idx_lo, idx_hi = 0, len(bin_centers) - 1

    range_lo = float(bin_centers[idx_lo])
    range_hi = float(bin_centers[idx_hi])
    print(
        f"Percentile range [{pct_lo}%–{pct_hi}%]: " f"[{range_lo:.3f}, {range_hi:.3f}]"
    )

    # Slice the histogram to the percentile-based range
    filtered_counts = hist_counts[idx_lo : idx_hi + 1]
    filtered_edges = bin_edges[idx_lo : idx_hi + 2]

    # Compute Multi-Otsu thresholds
    thresholds, fig = compute_otsu_thresholds(
        filtered_counts,
        classes=snakemake.params.otsu_k,
        bin_edges=filtered_edges,
        return_figure=True,
    )
    print(f"Batch Otsu thresholds: {[f'{t:.3f}' for t in thresholds]}")

    fig.savefig(snakemake.output.thresholds_png)

    # Save thresholds to TSV
    with open(snakemake.output.thresholds_tsv, "w") as fh:
        fh.write("threshold_index\tthreshold_value\n")
        for idx, val in enumerate(thresholds):
            fh.write(f"{idx}\t{val:.6f}\n")

    print(f"Saved thresholds to {snakemake.output.thresholds_tsv}")
    print(f"Saved histogram PNG to {snakemake.output.thresholds_png}")


if __name__ == "__main__":
    main()
