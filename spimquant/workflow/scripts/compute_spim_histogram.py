"""Compute intensity histogram from a raw SPIM dataset channel.

Reads the raw OME-Zarr at a given pyramid level using Dask and computes an
intensity histogram over a fixed range.  The resulting bin edges and counts
are saved as a compressed NumPy archive (.npz) for downstream aggregation
and threshold computation.

This is a Snakemake script; the ``snakemake`` object is automatically provided
when executed as part of a Snakemake workflow.
"""

import numpy as np
from dask_setup import get_dask_client
from zarrnii import ZarrNii


def main():
    level = snakemake.params.level
    stain = snakemake.wildcards.stain
    hist_bins = snakemake.params.hist_bins
    hist_range = snakemake.params.hist_range

    print(f"Computing histogram for stain={stain} at level={level}")
    print(f"  hist_bins={hist_bins}, hist_range={hist_range}")

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

    hist_counts = np.asarray(hist_counts, dtype=np.float64)
    bin_edges = np.asarray(bin_edges, dtype=np.float64)

    np.savez_compressed(
        snakemake.output.histogram,
        counts=hist_counts,
        bin_edges=bin_edges,
    )
    print(f"Saved histogram to {snakemake.output.histogram}")


if __name__ == "__main__":
    main()
