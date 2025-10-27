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
    from zarrnii.analysis import compute_otsu_thresholds

    # we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.corrected, **snakemake.params.zarrnii_kwargs
    )

    # first calculate histogram - using preset bins to avoid issues where bins are too large
    # because of high intensity outliers
    (hist_counts, bin_edges) = znimg.compute_histogram(bins=1000, range=(0, 1000))

    # get otsu thresholds (uses histogram)
    print("computing thresholds")
    (thresholds, fig) = compute_otsu_thresholds(
        hist_counts,
        classes=snakemake.params.otsu_k,
        bin_edges=bin_edges,
        return_figure=True,
    )
    print(f"  📈 thresholds: {[f'{t:.3f}' for t in thresholds]}")

    fig.savefig(snakemake.output.thresholds_png)

    print("thresholding image, saving as ome zarr")
    znimg_mask = znimg.segment_threshold(
        thresholds[snakemake.params.otsu_threshold_index]
    )

    # multiplying binary mask by 100 (so values are 0  and 100) to enable
    # field fraction calculation by subsequent local-mean downsampling
    znimg_mask = znimg_mask * 100

    # write to ome_zarr
    znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5)
