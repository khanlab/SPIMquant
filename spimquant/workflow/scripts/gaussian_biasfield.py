import dask
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii
from zarrnii.plugins import GaussianBiasFieldCorrection

hires_level = int(snakemake.wildcards.level)
proc_level = int(snakemake.params.proc_level)

unadjusted_downsample_factor = 2**proc_level
adjusted_downsample_factor = unadjusted_downsample_factor / (2**hires_level)


znimg = ZarrNii.from_ome_zarr(
    snakemake.input.spim,
    channel_labels=[snakemake.wildcards.stain],
    level=hires_level,
    downsample_near_isotropic=True,
    **snakemake.params.zarrnii_kwargs,
)

dask.config.set(scheduler="threads", num_workers=snakemake.threads)

print("compute bias field correction")
with ProgressBar():

    # Apply bias field correction
    znimg_corrected = znimg.apply_scaled_processing(
        GaussianBiasFieldCorrection(sigma=5.0),
        downsample_factor=adjusted_downsample_factor,
        upsampled_ome_zarr_path=snakemake.output.biasfield,
    )

    # write to ome_zarr
    znimg_corrected.to_ome_zarr(snakemake.output.corrected, max_layer=5)
