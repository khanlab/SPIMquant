"""Compute centroids from segmentation masks using ZarrNii.

This script reads a segmentation mask from an OME-Zarr file and computes the
centroids of labeled regions using ZarrNii's compute_centroids method. The
resulting centroids are saved as a NumPy array.

Example Dask cluster setup for advanced users:
"""

if __name__ == "__main__":
    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=int(snakemake.threads / 2),  # or 32, depending on workload
        threads_per_worker=2,  # isolate GIL
        memory_limit="auto",  # or tune to your RAM
        dashboard_address=":8788",
    )
    client = Client(cluster)
    print(cluster.dashboard_link)

    from zarrnii import ZarrNii
    import numpy as np
    import pandas as pd

    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.mask,
        level=0,  # input image is already downsampled to the wildcard level
        **snakemake.params.zarrnii_kwargs,
    )

    centroids = znimg.compute_centroids(output_path=snakemake.output.centroids_parquet)
