"""Aggregate intensity histograms across subjects.

Loads histogram ``.npz`` files (each containing ``counts`` and ``bin_edges``)
from multiple subjects and produces a single combined histogram by summing
counts.  All input histograms must share the same bin edges (i.e. they must
have been computed with the same ``hist_range`` and ``hist_bins`` parameters).

This is a Snakemake script; the ``snakemake`` object is automatically provided
when executed as part of a Snakemake workflow.
"""

import numpy as np


def main():
    histogram_files = snakemake.input.histograms

    print(f"Aggregating {len(histogram_files)} histogram(s)")

    total_counts = None
    ref_bin_edges = None

    for hist_file in histogram_files:
        data = np.load(hist_file)
        counts = data["counts"].astype(np.float64)
        bin_edges = data["bin_edges"]

        if ref_bin_edges is None:
            ref_bin_edges = bin_edges
            total_counts = counts.copy()
        else:
            if counts.shape == total_counts.shape and np.allclose(
                bin_edges, ref_bin_edges
            ):
                total_counts += counts
            else:
                # Rebin by treating bin centres as weighted samples
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
                rebinned, _ = np.histogram(
                    bin_centers,
                    bins=ref_bin_edges,
                    weights=counts,
                )
                total_counts += rebinned

        print(f"  {hist_file}: total voxels = {counts.sum():.0f}")

    print(f"Aggregated total voxels = {total_counts.sum():.0f}")

    np.savez_compressed(
        snakemake.output.histogram,
        counts=total_counts,
        bin_edges=ref_bin_edges,
    )
    print(f"Saved aggregated histogram to {snakemake.output.histogram}")


if __name__ == "__main__":
    main()
