from dask.diagnostics import ProgressBar

from zarrnii import ZarrNii


znimg = ZarrNii.from_ome_zarr(
    snakemake.input.spim,
    level=int(snakemake.wildcards.level),
    channel_labels=[snakemake.wildcards.stain],
    downsample_near_isotropic=True,
    **snakemake.params.zarrnii_kwargs,
)

with ProgressBar():
    znimg.to_nifti(snakemake.output.nii)
