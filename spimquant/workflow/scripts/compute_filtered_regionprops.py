"""Compute region properties from filtered segmentation masks using ZarrNii.
This script reads a segmentation mask from an OME-Zarr file, performs
connected components on chunks with overlap, applies filters based on
region properties, and outputs region properties on these filtered objects
"""

import tempfile
import zipfile
from contextlib import contextmanager

import numpy as np
import pandas as pd
from dask_setup import get_dask_client
from zarrnii import ZarrNii


@contextmanager
def get_zarr_path(mask_path):
    """Yield the path to the OME-Zarr store, extracting from zip if needed."""
    if mask_path.endswith(".ozx") or mask_path.endswith(".zip"):
        with tempfile.TemporaryDirectory(suffix=".ome.zarr") as temp_dir:
            print(f"Extracting zip archive to temporary directory: {temp_dir}")
            with zipfile.ZipFile(mask_path, "r") as zip_ref:
                zip_ref.extractall(temp_dir)
            yield temp_dir
    else:
        yield mask_path


if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        with get_zarr_path(snakemake.input.mask) as zarr_path:
            znimg = ZarrNii.from_file(zarr_path, level=0)

            znimg.compute_region_properties(
                output_path=snakemake.output.regionprops_parquet,
                region_filters=snakemake.params.region_filters,
                output_properties=snakemake.params.output_properties,
            )
            voxel_volume_mm3 = float(np.prod([float(z) for z in znimg.get_zooms()]))

    # physical object volume, so sizes compare across scans and resolutions
    nvox = snakemake.params.output_properties.get("area", "area")
    df = pd.read_parquet(snakemake.output.regionprops_parquet)
    df["volume"] = df[nvox].astype("float64") * voxel_volume_mm3
    df.to_parquet(snakemake.output.regionprops_parquet, index=False)
