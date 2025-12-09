"""
Merge segmentation statistic TSV tables, drop coordinate columns, and aggregate.

This script loads multiple TSV files containing segmentation statistics at
instance and aggregate levels. All tables are merged on ('index', 'name')
using an outer join. Columns ending with '_x', '_y', or '_z' are dropped.

After merging, rows sharing the same ('index', 'name') are collapsed.
Aggregation rules:
  - 'nvoxels' is summed across matching entries.
  - All other numeric columns are averaged.
  - Non-numeric columns (other than index/name) are preserved by first-occurrence.

Density is recomputed as count/volume if both are available.
Finally, the merged table is written to output.
"""

import pandas as pd
from functools import reduce

# read each TSV into a list of dataframes
dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input]

# merge all on index + name
merged = reduce(
    lambda left, right: pd.merge(left, right, on=("index", "name"), how="outer"), dfs
)

# ---- drop unwanted coordinate columns ----
drop_cols = [c for c in merged.columns if c.endswith(("_x", "_y", "_z"))]
merged = merged.drop(columns=drop_cols)

# fill all columns except 'volume'
cols_to_fill = [c for c in merged.columns if c != "volume"]
merged[cols_to_fill] = merged[cols_to_fill].fillna(0)

# ---- aggregation rules ----
agg_dict = {}
for col in merged.columns:
    if col in ["index", "name"]:
        continue
    elif col == "nvoxels":
        agg_dict[col] = "sum"
    else:
        # mean for all other numeric columns,
        # 'first' for non-numeric to avoid errors
        agg_dict[col] = (
            "mean" if pd.api.types.is_numeric_dtype(merged[col]) else "first"
        )

merged = merged.groupby(["index", "name"], as_index=False).agg(agg_dict)

# recompute density if available
if "count" in merged.columns and "volume" in merged.columns:
    merged["density"] = merged["count"] / merged["volume"]

if hasattr(snakemake.params, "columns_to_drop"):
    merged = merged.drop(columns=snakemake.params.columns_to_drop)

merged.to_csv(snakemake.output.tsv, sep="\t", index=False)
