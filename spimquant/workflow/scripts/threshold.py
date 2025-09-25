from zarrnii import ZarrNii

znimg_hires = ZarrNii.from_ome_zarr(snakemake.input.corrected)

znimg_hires.segment_threshold(snakemake.params.threshold).to_ome_zarr(
    snakemake.output.mask, max_layer=5
)
