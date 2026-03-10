if __name__ == "__main__":

    from dask_setup import get_dask_client
    from zarrnii import ZarrNii
    from zarrnii.plugins import N4BiasFieldCorrection

    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

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

        print("compute bias field correction")

        adjusted_chunk = int(320 / (2**adjusted_downsample_factor))

        # Apply bias field correction
        znimg_corrected = znimg.apply_scaled_processing(
            N4BiasFieldCorrection(shrink_factor=snakemake.params.shrink_factor),
            downsample_factor=adjusted_downsample_factor,
            upsampled_ome_zarr_path=snakemake.output.biasfield,
        )

        # write to ome_zarr
        znimg_corrected.to_ome_zarr(snakemake.output.corrected, max_layer=5)
