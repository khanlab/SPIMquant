"""Aggregate transformed regionprops across subjects.

This script combines transformed regionprops from multiple subjects into a single
table, adding a 'subject' column to identify the source of each row.
"""

import pandas as pd

# Load all the aggregated regionprops parquet files (already aggregated across stains)
dfs = []
for subject, parquet_file in zip(snakemake.params.subjects, snakemake.input.regionprops_parquets):
    df = pd.read_parquet(parquet_file)
    df['subject'] = subject
    dfs.append(df)

# Concatenate all dataframes
df_combined = pd.concat(dfs, ignore_index=True)

# Save the aggregated output
df_combined.to_parquet(snakemake.output.regionprops_aggregated_parquet, index=False)
