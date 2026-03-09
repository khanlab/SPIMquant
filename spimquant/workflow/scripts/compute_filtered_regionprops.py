"""Compute region properties from filtered segmentation masks using ZarrNii.

This script reads a segmentation mask from an OME-Zarr file, performs
connected components on chunks with overlap, applys filters based on 
region properties, and outputs region properties on these filtered objects
"""

from dask_setup import get_dask_client
from zarrnii import ZarrNii

with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

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
