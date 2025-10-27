if __name__ == "__main__":

    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=int(snakemake.threads / 2),  # or 32, depending on workload
        threads_per_worker=2,  # isolate GIL
        memory_limit="auto",  # or tune to your RAM
        dashboard_address=":8788",
    )
    client = Client(cluster)
    print(cluster.dashboard_link)

    from zarrnii import ZarrNii
    from zarrnii.plugins import N4BiasFieldCorrection

    hires_level = int(snakemake.wildcards.level)
    ds_level = int(snakemake.wildcards.dslevel)

    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.spim,
        channel_labels=[snakemake.wildcards.stain],
        level=hires_level,
        **snakemake.params.zarrnii_kwargs,
    )

    print("compute bias field correction")

    # Apply bias field correction
    znimg_corrected = znimg.apply_scaled_processing(
        N4BiasFieldCorrection(sigma=5.0),
        downsample_factor=2**ds_level,
        upsampled_ome_zarr_path=snakemake.output.biasfield,
    )

    # write to ome_zarr
    znimg_corrected.to_ome_zarr(snakemake.output.corrected, max_layer=1)
