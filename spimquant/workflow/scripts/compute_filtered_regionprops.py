"""Compute region properties from filtered segmentation masks using ZarrNii.

This script reads a segmentation mask from an OME-Zarr file, performs
connected components on chunks with overlap, applys filters based on 
region properties, and outputs region properties on these filtered objects
"""

if __name__ == "__main__":
    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=int(snakemake.threads / 4),  # or 32, depending on workload
        threads_per_worker=4,  # isolate GIL
        memory_limit="auto",  # or tune to your RAM
        dashboard_address=":8788",
    )
    client = Client(cluster)
    print(cluster.dashboard_link)

    try:
        from zarrnii import ZarrNii
        import numpy as np
        import pandas as pd

        znimg = ZarrNii.from_ome_zarr(
            snakemake.input.mask,
            level=0,  # input image is already downsampled to the wildcard level
            **snakemake.params.zarrnii_kwargs,
        )

        znimg.compute_region_properties(
            output_path=snakemake.output.regionprops_parquet,
            region_filters=snakemake.params.region_filters,
            output_properties=snakemake.params.output_properties,
        )
    finally:
        client.close()
        cluster.close()
