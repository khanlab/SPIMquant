import numpy as np

from dask_setup import get_dask_client
from zarrnii import ZarrNii
from zarrnii.analysis import compute_otsu_thresholds
import matplotlib

matplotlib.use("agg")

if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

        zarrnii_kwargs = snakemake.params.zarrnii_kwargs
        pct_lo, pct_hi = snakemake.params.hist_percentile_range
        bin_width = snakemake.params.hist_bin_width

        # load a downsampled version to estimate the percentile-based range
        print("estimating intensity range from downsampled image...")
        znimg_ds = None
        for ds_level in [5, 4, 3, 2, 1]:
            try:
                candidate = ZarrNii.from_ome_zarr(
                    snakemake.input.corrected, level=ds_level, **zarrnii_kwargs
                )
                znimg_ds = candidate
                break
            except Exception:
                pass

        if znimg_ds is None:
            znimg_ds = ZarrNii.from_ome_zarr(
                snakemake.input.corrected, **zarrnii_kwargs
            )

        data_ds = znimg_ds.data.compute().ravel().astype(np.float32)
        range_lo = float(np.percentile(data_ds, pct_lo))
        range_hi = float(np.percentile(data_ds, pct_hi))
        print(
            f"  📊 percentile range [{pct_lo}%, {pct_hi}%]: [{range_lo:.3f}, {range_hi:.3f}]"
        )

        # compute number of bins from bin width
        n_bins = max(2, int(np.ceil((range_hi - range_lo) / bin_width)))
        print(f"  📊 bins: {n_bins} (bin width: {bin_width})")

        # we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
        znimg = ZarrNii.from_ome_zarr(snakemake.input.corrected, **zarrnii_kwargs)

        # calculate histogram using percentile-based range and bin-width-derived bin count
        (hist_counts, bin_edges) = znimg.compute_histogram(
            bins=n_bins, range=[range_lo, range_hi]
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
