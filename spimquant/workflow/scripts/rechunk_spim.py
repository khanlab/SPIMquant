from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from zarrnii import ZarrNii

if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        chunk_size_zyx = tuple(int(x) for x in snakemake.params.chunk_size)
        znimg = ZarrNii.from_file(
            snakemake.input.spim,
            level=int(snakemake.wildcards.level),
            downsample_near_isotropic=True,
            **snakemake.params.zarrnii_kwargs,
        )

        if znimg.darr.ndim < 3:
            raise ValueError(
                f"Expected at least 3 dimensions for rechunking, got {znimg.darr.ndim}."
            )

        leading_chunks = tuple(znimg.darr.chunksize[:-3])
        target_chunks = leading_chunks + chunk_size_zyx
        znimg.darr = znimg.darr.rechunk(target_chunks)

        with ProgressBar():
            znimg.to_ome_zarr(
                snakemake.output.rechunked,
                match_scale_factors_from=snakemake.input.spim,
            )
