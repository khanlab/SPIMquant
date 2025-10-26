if __name__ == "__main__":

    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=int(snakemake.threads/2),             # or 32, depending on workload
        threads_per_worker=2,     # isolate GIL
        memory_limit="auto",       # or tune to your RAM
        dashboard_address=':8788',
    )
    client = Client(cluster)
    print(cluster.dashboard_link)


    from zarrnii import ZarrNii

    # we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.corrected, **snakemake.params.zarrnii_kwargs
    )


    # get otsu thresholds (uses histogram)
    print("computing thresholds")
    thresholds = znimg.compute_otsu_thresholds(classes=snakemake.params.otsu_k)
    print(f"  📈 thresholds: {[f'{t:.3f}' for t in thresholds]}")

    print("thresholding image, saving as ome zarr")
    znimg_mask = znimg.segment_threshold(
        thresholds[snakemake.params.otsu_threshold_index]
    )

    # write to ome_zarr
    znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5)
