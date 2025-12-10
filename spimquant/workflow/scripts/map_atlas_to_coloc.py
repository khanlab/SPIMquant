"""Map computed coloc objects to atlas regions and generate labeled statistics.
"""

from zarrnii import ZarrNiiAtlas
import pandas as pd

coloc = pd.read_parquet(snakemake.input.coloc_parquet).to_dict(orient="list")
atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

df_coloc, df_counts = atlas.label_region_properties(
    coloc,
    coord_column_names=snakemake.params.coord_column_names,
    include_names=True,
)

df_coloc.to_csv(snakemake.output.coloc_tsv, sep="\t", index=False)
df_counts.to_csv(snakemake.output.counts_tsv, sep="\t", index=False)
