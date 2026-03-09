from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from zarrnii import ZarrNii

with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.spim,
        level=int(snakemake.wildcards.level),
        channel_labels=[snakemake.wildcards.stain],
        downsample_near_isotropic=True,
        **snakemake.params.zarrnii_kwargs,
    )

    with ProgressBar():
        znimg.to_nifti(snakemake.output.nii)
