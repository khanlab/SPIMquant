"""Merge multiple segmentation statistics TSV files.

This script merges multiple TSV files containing segmentation statistics
(e.g., counts, volumes) for different regions. The files are merged on the
'index' and 'name' columns using an outer join. After merging, density metrics
are computed as count/volume when both columns are present.
"""

import pandas as pd
from functools import reduce

# read each TSV into a list of dataframes
dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input]

# merge all dataframes on 'index'

merged = reduce(
    lambda left, right: pd.merge(left, right, on=("index", "name"), how="outer"), dfs
)

# fill all columns except volume with 0
cols_to_fill = [c for c in merged.columns if c != "volume"]
merged[cols_to_fill] = merged[cols_to_fill].fillna(0)

# compute density only if columns exist
if "count" in merged.columns and "volume" in merged.columns:
    merged["density"] = merged["count"] / merged["volume"]

if hasattr(snakemake.params, "columns_to_drop"):
    merged = merged.drop(columns=snakemake.params.columns_to_drop)

merged.to_csv(snakemake.output.tsv, sep="\t", index=False)
