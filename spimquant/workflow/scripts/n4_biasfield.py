from dask.diagnostics import ProgressBar

from zarrnii.plugins import N4BiasFieldCorrection
from zarrnii import ZarrNii

hires_level = int(snakemake.wildcards.level)
ds_level = int(snakemake.wildcards.dslevel)

znimg = ZarrNii.from_ome_zarr(
    snakemake.input.spim,
    channel_labels=[snakemake.wildcards.stain],
    level=hires_level,
    downsample_near_isotropic=True,
)


print("compute bias field correction")
with ProgressBar():

    # Apply bias field correction
    znimg_corrected = znimg.apply_scaled_processing(
        N4BiasFieldCorrection(sigma=5.0),
        downsample_factor=2**ds_level,
        upsampled_ome_zarr_path=snakemake.output.biasfield,
    )

    # write to ome_zarr
    znimg_corrected.to_ome_zarr(snakemake.output.corrected, max_layer=5)
