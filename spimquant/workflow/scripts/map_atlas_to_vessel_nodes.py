"""Map vessel node coordinates to atlas regions and generate labeled statistics.

This script takes vessel node coordinates (from a nodes parquet file) and maps
them to atlas regions using a ZarrNiiAtlas object, generating two outputs:
1. Annotated nodes parquet with atlas region labels assigned to each node
2. Count statistics showing the number of nodes per atlas region
"""

import pandas as pd
from zarrnii import ZarrNiiAtlas

nodes = pd.read_parquet(snakemake.input.nodes_parquet).to_dict(orient="list")
atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

df_nodes, df_counts = atlas.label_region_properties(
    nodes,
    coord_column_names=snakemake.params.coord_column_names,
    include_names=True,
)

df_nodes.to_parquet(snakemake.output.nodes_parquet, index=False)
df_counts.to_csv(snakemake.output.counts_tsv, sep="\t", index=False)
