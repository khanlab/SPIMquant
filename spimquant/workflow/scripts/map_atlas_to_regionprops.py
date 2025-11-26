"""Map computed regionprops to atlas regions and generate labeled statistics.

This script takes computed regionprops from segmentation and maps them to atlas
regions, generating two outputs:
1. Labeled objects with their regionprops and corresponding atlas label assignments
2. Count statistics showing the number of objects per atlas region
"""

from zarrnii import ZarrNiiAtlas
import pandas as pd

regionprops = pd.read_parquet(snakemake.input.regionprops_parquet).to_dict(
    orient="list"
)

atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

df_regionprops, df_counts = atlas.label_region_properties(
    regionprops,
    coord_column_names=snakemake.params.coord_column_names,
    include_names=True,
)

df_regionprops.to_csv(snakemake.output.regionprops_tsv, sep="\t", index=False)
df_counts.to_csv(snakemake.output.counts_tsv, sep="\t", index=False)
