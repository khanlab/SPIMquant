if __name__ == "__main__":
    """
    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=2, #int(snakemake.threads / 2),  # or 32, depending on workload
        threads_per_worker=2,  # isolate GIL
        memory_limit="12GB", #"auto",  # or tune to your RAM
        dashboard_address=":8788",
    )
    client = Client(cluster)
    print(cluster.dashboard_link)
    """
    from dask.diagnostics import ProgressBar

    from zarrnii import ZarrNii
    import numpy as np

    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.mask,
        level=0,
        **snakemake.params.zarrnii_kwargs,
    )
    with ProgressBar():
        centroids = znimg.compute_centroids()

    np.save(snakemake.output.centroids_npy, centroids)
