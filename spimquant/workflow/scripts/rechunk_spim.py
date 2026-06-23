"""Rechunk SPIM OME-Zarr to work directory for optimized processing.
"""

from zarrnii import ZarrNii
from dask_setup import get_dask_client

if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        znimg = ZarrNii.from_file(
            snakemake.input.spim,
            chunks=snakemake.params.chunks,
            **snakemake.params.zarrnii_kwargs,
        )
        znimg.to_ome_zarr(snakemake.output.spim)
