"""Compute and save an intensity histogram for a single subject.

Used as the first step of group-level Otsu thresholding.  Each subject
independently computes its histogram from the bias-field corrected image;
the resulting per-subject NPZ files are later aggregated by the
``group_otsu`` rule to derive a single set of thresholds shared across the
whole cohort.

This is a Snakemake script; the ``snakemake`` object is automatically
provided when executed as part of a Snakemake workflow.
"""

import numpy as np

from dask_setup import get_dask_client
from zarrnii import ZarrNii

if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

        zarrnii_kwargs = snakemake.params.zarrnii_kwargs
        pct_lo, pct_hi = snakemake.params.hist_percentile_range
        bin_width = snakemake.params.hist_bin_width

        # load a downsampled version to estimate the percentile-based range
        print("estimating intensity range from downsampled image...")
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
        print(
            f"  📊 percentile range [{pct_lo}%, {pct_hi}%]: [{range_lo:.3f}, {range_hi:.3f}]"
        )

        # compute number of bins from bin width
        n_bins = max(2, int(np.ceil((range_hi - range_lo) / bin_width)))
        print(f"  📊 bins: {n_bins} (bin width: {bin_width})")

        # compute full-resolution histogram
        znimg = ZarrNii.from_ome_zarr(snakemake.input.corrected, **zarrnii_kwargs)
        (hist_counts, bin_edges) = znimg.compute_histogram(
            bins=n_bins, range=[range_lo, range_hi]
        )

        print(f"saving histogram to {snakemake.output.histogram_npz}")
        np.savez(
            snakemake.output.histogram_npz,
            hist_counts=hist_counts,
            bin_edges=bin_edges,
        )
