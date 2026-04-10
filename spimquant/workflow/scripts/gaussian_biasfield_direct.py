import dask.array as da
import numpy as np
from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from zarrnii import ZarrNii
from zarrnii.plugins import GaussianBiasFieldCorrection

if __name__ == "__main__":
    hires_level = int(snakemake.wildcards.level)

    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

        znimg = ZarrNii.from_ome_zarr(
            snakemake.input.spim,
            channel_labels=[snakemake.wildcards.stain],
            level=hires_level,
            downsample_near_isotropic=True,
            **snakemake.params.zarrnii_kwargs,
        )

        print("compute bias field correction (direct, no downsample/upsample)")
        plugin = GaussianBiasFieldCorrection(sigma=5.0)

        with ProgressBar():
            # Compute bias field directly at the current level without downsampling
            array = znimg.data.compute()
            bias_field = plugin.lowres_func(array)

            # Write bias field to ome_zarr
            bias_znimg = znimg.copy()
            bias_znimg.data = da.from_array(bias_field, chunks=znimg.data.chunks)
            bias_znimg.to_ome_zarr(snakemake.output.biasfield, max_layer=5)

            # Apply correction directly (bias field at same resolution as input)
            corrected_znimg = znimg.copy()
            corrected_znimg.data = plugin.highres_func(znimg.data, bias_znimg.data)

            # write to ome_zarr
            corrected_znimg.to_ome_zarr(snakemake.output.corrected, max_layer=5)
