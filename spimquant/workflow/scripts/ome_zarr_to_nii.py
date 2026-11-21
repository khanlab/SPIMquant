import remote_zarr  # noqa: F401  OME-Zarr outputs as pointer files (--zarr_store_root)
from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from zarrnii import ZarrNii

if __name__ == "__main__":
    with get_dask_client("threads", snakemake.threads):
        znimg = ZarrNii.from_file(
            snakemake.params.spim,
            level=int(snakemake.wildcards.level),
            channel_labels=[snakemake.wildcards.stain],
            downsample_near_isotropic=True,
            chunks=None,
            **snakemake.params.zarrnii_kwargs,
        )

        with ProgressBar():
            znimg.to_nifti(snakemake.output.nii)
