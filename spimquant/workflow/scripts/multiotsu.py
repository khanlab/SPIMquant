from dask.diagnostics import ProgressBar

from zarrnii import ZarrNii

# we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
znimg = ZarrNii.from_ome_zarr(snakemake.input.corrected)


# get otsu thresholds (uses histogram)
print("computing thresholds")
with ProgressBar():
    thresholds = znimg.compute_otsu_thresholds(classes=snakemake.params.otsu_k)
print(f"  📈 thresholds: {[f'{t:.3f}' for t in thresholds]}")

print("thresholding image, saving as ome zarr")
with ProgressBar():
    znimg_mask = znimg.segment_threshold(
        thresholds[snakemake.params.otsu_threshold_index]
    )

    # write to ome_zarr
    znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5)
