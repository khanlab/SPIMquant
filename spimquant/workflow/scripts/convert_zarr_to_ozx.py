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



        # import ngff_zarr as nz
        from zarrnii import ZarrNii

        # Read from directory store
        # multiscales = nz.from_ngff_zarr(snakemake.input.zarr)

        # print(multiscales)
        # Write as .ozx file
        # nz.to_ngff_zarr(snakemake.output.ozx, multiscales, version='0.5')

        # note, this recomputes the multiscales
        ZarrNii.from_ome_zarr(snakemake.input.zarr).to_ome_zarr(snakemake.output.ozx)


    finally:
        client.close()
        cluster.close()
