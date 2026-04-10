import dask.array as da
from dask_setup import get_dask_client
from zarrnii import ZarrNii
from zarrnii.plugins import N4BiasFieldCorrection

if __name__ == "__main__":

    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

        hires_level = int(snakemake.wildcards.level)

        znimg = ZarrNii.from_ome_zarr(
            snakemake.input.spim,
            channel_labels=[snakemake.wildcards.stain],
            level=hires_level,
            downsample_near_isotropic=True,
            **snakemake.params.zarrnii_kwargs,
        )

        print("compute N4 bias field correction (direct, no downsample/upsample)")
        plugin = N4BiasFieldCorrection(shrink_factor=snakemake.params.shrink_factor)

        # Compute bias field directly at the current level without downsampling
        array = znimg.data.compute()
        bias_field = plugin.lowres_func(array)

        # Apply correction directly (bias field at same resolution as input)
        bias_field_dask = da.from_array(bias_field, chunks=znimg.data.chunks)
        corrected_znimg = znimg.copy()
        corrected_znimg.data = plugin.highres_func(znimg.data, bias_field_dask)

        # write to ome_zarr
        corrected_znimg.to_ome_zarr(snakemake.output.corrected, max_layer=5)
