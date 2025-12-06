"""Aggregate transformed regionprops across stains.

This script combines transformed regionprops from multiple stains into a single
table, adding a 'stain' column to identify the source of each row.
"""

import pandas as pd

# Load all the transformed regionprops parquet files
dfs = []
for stain, parquet_file in zip(snakemake.params.stains, snakemake.input.regionprops_parquets):
    df = pd.read_parquet(parquet_file)
    df['stain'] = stain
    dfs.append(df)

# Concatenate all dataframes
df_combined = pd.concat(dfs, ignore_index=True)

# Save the aggregated output
df_combined.to_parquet(snakemake.output.regionprops_aggregated_parquet, index=False)
