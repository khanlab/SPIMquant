if __name__ == "__main__":

    from dask_setup import get_dask_client
    from zarrnii import ZarrNii
    from zarrnii.plugins import N4BiasFieldApply

    is_imaris = str(snakemake.input.spim).lower().endswith(".ims")

    with get_dask_client(
        snakemake.config["dask_scheduler"],
        snakemake.threads,
        threads_per_worker=16 if is_imaris else 2,
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
            chunks=(256, 256, 256) if is_imaris else None,
            **snakemake.params.zarrnii_kwargs,
        )
        znimg_lowres = ZarrNii.from_nifti(snakemake.input.biasfield, axes_order="ZYX")

        print("compute bias field correction")
        scaled_proc_kwargs = {"lowres_znimg": znimg_lowres, "method": "map_blocks"}

        # scaled_proc_kwargs controls how apply_scaled_processing is performed

        # Apply bias field correction
        znimg_corrected = znimg.apply_scaled_processing(
            N4BiasFieldApply(log_space=False),
            **scaled_proc_kwargs,
        )

        # write to ome_zarr
        znimg_corrected.to_ome_zarr(
            snakemake.output.corrected,
            match_scale_factors_from=snakemake.input.spim,
            **snakemake.config["zarrnii_out_kwargs"],
        )
