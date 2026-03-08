from dask_setup import get_dask_client
from zarrnii import ZarrNii
from zarrnii.analysis import compute_otsu_thresholds
import matplotlib

matplotlib.use("agg")

with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

    # we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.corrected, **snakemake.params.zarrnii_kwargs
    )

    # first calculate histogram - using preset bins to avoid issues where bins are too large
    # because of high intensity outliers
    (hist_counts, bin_edges) = znimg.compute_histogram(
        bins=snakemake.params.hist_bins, range=snakemake.params.hist_range
    )

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
