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

    try:

        from zarrnii import ZarrNii

        # read to get shape
        znimg = ZarrNii.from_ome_zarr(
            snakemake.input.spim,
            **snakemake.params.zarrnii_kwargs,
        )
        shape = znimg.data.shape

        # read with re-chunk
        znimg = ZarrNii.from_ome_zarr(
            snakemake.input.spim,
            chunks=(1, 1, shape[-2], shape[-1]),
            rechunk=True,
            **snakemake.params.zarrnii_kwargs,
        )

        znimg_destriped = znimg.destripe(**snakemake.params.destripe_kwargs)

        # write to ome_zarr
        znimg_destriped.to_ome_zarr(snakemake.output.destriped, max_layer=0)

    finally:
        client.close()
        cluster.close()
