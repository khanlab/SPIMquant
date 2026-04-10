"""Aggregate per-subject histograms and compute group-level Otsu thresholds.

Loads the intensity histogram NPZ files produced by ``compute_subject_histogram``
for every subject, merges them onto a common intensity grid, and applies
multi-level Otsu thresholding to the aggregate histogram.  The resulting
thresholds are saved as a JSON file (for downstream use by ``multiotsu_group``)
and as a PNG figure for visual inspection.

This is a Snakemake script; the ``snakemake`` object is automatically
provided when executed as part of a Snakemake workflow.
"""

import json

import matplotlib

matplotlib.use("agg")
import numpy as np

from zarrnii.analysis import compute_otsu_thresholds

if __name__ == "__main__":
    # Load all per-subject histograms
    histograms = []
    for path in snakemake.input.histogram_npz:
        data = np.load(path)
        histograms.append((data["hist_counts"], data["bin_edges"]))

    print(f"Loaded {len(histograms)} subject histograms")

    # Find the common intensity range spanning all subjects
    overall_lo = min(float(be[0]) for _, be in histograms)
    overall_hi = max(float(be[-1]) for _, be in histograms)

    # Use the first histogram's bin width as reference for the common grid
    ref_bin_edges = histograms[0][1]
    bin_width = float(ref_bin_edges[1] - ref_bin_edges[0])

    n_bins = max(2, int(np.ceil((overall_hi - overall_lo) / bin_width)))
    common_bin_edges = np.linspace(overall_lo, overall_hi, n_bins + 1)

    print(
        f"  📊 common range: [{overall_lo:.3f}, {overall_hi:.3f}], bins: {n_bins}"
    )

    # Aggregate all subject histograms onto the common grid
    aggregate_counts = np.zeros(n_bins, dtype=np.float64)
    for hist_counts, bin_edges in histograms:
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        # Map each original bin center to the nearest common bin
        common_indices = np.searchsorted(common_bin_edges[1:], bin_centers)
        common_indices = np.clip(common_indices, 0, n_bins - 1)
        np.add.at(aggregate_counts, common_indices, hist_counts)

    # Apply multi-level Otsu thresholding to the aggregated histogram
    print("computing group Otsu thresholds")
    (thresholds, fig) = compute_otsu_thresholds(
        aggregate_counts,
        classes=snakemake.params.otsu_k,
        bin_edges=common_bin_edges,
        return_figure=True,
    )
    print(f"  📈 group thresholds: {[f'{t:.3f}' for t in thresholds]}")

    otsu_threshold_index = snakemake.params.otsu_threshold_index
    selected_threshold = float(thresholds[otsu_threshold_index])
    print(
        f"  ✅ selected threshold (index {otsu_threshold_index}): {selected_threshold:.3f}"
    )

    fig.savefig(snakemake.output.thresholds_png)

    # Save thresholds as JSON for use by the per-subject segmentation rule
    result = {
        "thresholds": thresholds.tolist(),
        "otsu_threshold_index": otsu_threshold_index,
        "selected_threshold": selected_threshold,
        "n_subjects": len(histograms),
    }
    with open(snakemake.output.thresholds_json, "w") as f:
        json.dump(result, f, indent=2)

    print(f"saved group thresholds to {snakemake.output.thresholds_json}")
