"""Apply a precomputed group-level Otsu threshold to a single subject.

Reads the group threshold JSON produced by ``group_otsu``, applies it to
the bias-field corrected image for the current subject, and writes the
binary mask.  A per-subject PNG is also saved that overlays the group
threshold on this subject's own intensity histogram, which is useful for
visual quality control.

This is a Snakemake script; the ``snakemake`` object is automatically
provided when executed as part of a Snakemake workflow.
"""

import json

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np

from dask_setup import get_dask_client
from zarrnii import ZarrNii

if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

        # Load group threshold from JSON
        with open(snakemake.input.thresholds_json) as f:
            group_data = json.load(f)

        all_thresholds = group_data["thresholds"]
        otsu_threshold_index = group_data["otsu_threshold_index"]
        selected_threshold = group_data["selected_threshold"]

        print(f"  📈 group thresholds: {[f'{t:.3f}' for t in all_thresholds]}")
        print(f"  ✅ applying threshold: {selected_threshold:.3f}")

        zarrnii_kwargs = snakemake.params.zarrnii_kwargs
        pct_lo, pct_hi = snakemake.params.hist_percentile_range
        bin_width = snakemake.params.hist_bin_width

        # Load a downsampled version to estimate the percentile-based range for
        # the per-subject histogram visualisation
        print("estimating intensity range from downsampled image for QC figure...")
        znimg_ds = None
        for ds_level in [5, 4, 3, 2, 1]:
            try:
                candidate = ZarrNii.from_ome_zarr(
                    snakemake.input.corrected, level=ds_level, **zarrnii_kwargs
                )
                znimg_ds = candidate
                break
            except Exception:
                pass

        if znimg_ds is None:
            znimg_ds = ZarrNii.from_ome_zarr(
                snakemake.input.corrected, **zarrnii_kwargs
            )

        data_ds = znimg_ds.data.compute().ravel().astype(np.float32)
        range_lo = float(np.percentile(data_ds, pct_lo))
        range_hi = float(np.percentile(data_ds, pct_hi))

        n_bins = max(2, int(np.ceil((range_hi - range_lo) / bin_width)))

        # Compute full-resolution histogram for this subject
        znimg = ZarrNii.from_ome_zarr(snakemake.input.corrected, **zarrnii_kwargs)
        (hist_counts, bin_edges) = znimg.compute_histogram(
            bins=n_bins, range=[range_lo, range_hi]
        )

        # Generate per-subject visualisation: subject histogram + group threshold
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.bar(
            bin_centers,
            hist_counts,
            width=bin_width,
            color="steelblue",
            alpha=0.7,
            label="Subject histogram",
        )
        for j, t in enumerate(all_thresholds):
            linestyle = "-" if j == otsu_threshold_index else "--"
            label = f"Group threshold[{j}]={t:.3f}"
            if j == otsu_threshold_index:
                label += " (selected)"
            ax.axvline(t, color="red", linestyle=linestyle, label=label)
        ax.set_xlabel("Intensity")
        ax.set_ylabel("Count")
        ax.set_title(
            f"Subject histogram with group Otsu threshold "
            f"(selected: {selected_threshold:.3f})"
        )
        ax.legend()
        fig.tight_layout()
        fig.savefig(snakemake.output.thresholds_png)
        plt.close(fig)

        # Apply the group threshold to create the binary mask
        print("thresholding image with group threshold, saving as ome zarr")
        znimg_mask = znimg.segment_threshold(selected_threshold)

        # Multiply by 100 (values 0 and 100) to enable field-fraction
        # calculation by subsequent local-mean downsampling
        znimg_mask = znimg_mask * 100

        znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5)
