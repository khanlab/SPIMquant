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
    from zarrnii.plugins import SegmentationCleaner

    hires_level = int(snakemake.wildcards.level)
    ds_level = int(snakemake.wildcards.dslevel)

    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.mask,
        level=0,
        **snakemake.params.zarrnii_kwargs,
    )

    # perform cleaning of artifactual positives by
    # removing objects with low extent (extent is ratio of num voxels to bounding box)

    # the downsample_factor we use should be proportional to the segmentation level
    #   e.g. if segmentation level is 3, then we have already downsampled by 2^3, so
    #   the downsample factor should be divided by that..
    unadjusted_downsample_factor = 2**ds_level

    adjusted_downsample_factor = unadjusted_downsample_factor / (2**hires_level)

    znimg_cleaned = znimg.apply_scaled_processing(
        SegmentationCleaner(max_extent=snakemake.params.max_extent),
        downsample_factor=adjusted_downsample_factor,
        upsampled_ome_zarr_path=snakemake.output.exclude_mask,
    )

    # write to final ome_zarr
    znimg_cleaned.to_ome_zarr(snakemake.output.cleaned_mask, max_layer=5)
