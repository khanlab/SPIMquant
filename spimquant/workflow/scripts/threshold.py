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

    znimg_hires = ZarrNii.from_ome_zarr(
        snakemake.input.corrected, **snakemake.params.zarrnii_kwargs
    )

    print("thresholding image, saving as ome zarr")
    znimg_mask = znimg.segment_threshold(snakemake.params.threshold)

    # multiplying binary mask by 100 (so values are 0  and 100) to enable
    # field fraction calculation by subsequent local-mean downsampling
    znimg_mask = znimg_mask * 100

    # write to ome_zarr
    znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5)
