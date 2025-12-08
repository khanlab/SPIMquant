import dask
from dask.diagnostics import ProgressBar

from zarrnii import ZarrNii

dask.config.set(scheduler="threads", num_workers=snakemake.threads)

znimg = ZarrNii.from_ome_zarr(
    snakemake.input.spim,
    level=int(snakemake.wildcards.level),
    channel_labels=[snakemake.wildcards.stain],
    downsample_near_isotropic=True,
    **snakemake.params.zarrnii_kwargs,
)

with ProgressBar():
    znimg.to_nifti(snakemake.output.nii)
