if __name__ == "__main__":

    from dask_setup import get_dask_client
    from zarrnii import ZarrNii
    from zarrnii.plugins import N4BiasFieldCorrection

    with get_dask_client(
        snakemake.config["dask_scheduler"], snakemake.threads, threads_per_worker=16
    ):

        hires_level = int(snakemake.wildcards.level)
        proc_level = int(snakemake.params.proc_level)

        unadjusted_downsample_factor = 2**proc_level
        adjusted_downsample_factor = unadjusted_downsample_factor / (2**hires_level)

        znimg = ZarrNii.from_file(
            snakemake.input.spim,
            channel_labels=[snakemake.wildcards.stain],
            level=hires_level,
            downsample_near_isotropic=True,
            chunks=(256, 256, 256),
            **snakemake.params.zarrnii_kwargs,
        )
        znimg_lowres = ZarrNii.from_file(
            snakemake.input.spim,
            channel_labels=[snakemake.wildcards.stain],
            level=proc_level,
            downsample_near_isotropic=True,
            **snakemake.params.zarrnii_kwargs,
        )

        print("compute bias field correction")

        # Apply bias field correction
        znimg_corrected = znimg.apply_scaled_processing(
            N4BiasFieldCorrection(shrink_factor=snakemake.params.shrink_factor),
            lowres_znimg=znimg_lowres,
            method="map_blocks",
        )

        # write to ome_zarr
        znimg_corrected.to_ome_zarr(
            snakemake.output.corrected, match_scale_factors_from=snakemake.input.spim
        )
