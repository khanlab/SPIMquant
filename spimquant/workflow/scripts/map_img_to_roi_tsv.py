from zarrnii import ZarrNii, ZarrNiiAtlas

atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)
img = ZarrNii.from_file(snakemake.input.img)

dseg_df = atlas.aggregate_image_by_regions(
    img, aggregation_func="mean", column_name=snakemake.wildcards.suffix
)

dseg_df.to_csv(snakemake.output.tsv, sep="\t", index=False)
