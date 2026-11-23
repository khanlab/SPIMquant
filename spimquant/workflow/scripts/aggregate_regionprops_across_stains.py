"""Aggregate the template-space objects of several stains into one table.

Each input already carries ``stain``, ``seg_method`` and ``object_method``
columns (added by transform_regionprops_to_template.py), so this only
concatenates them.
"""

import pandas as pd

dfs = [pd.read_parquet(f) for f in snakemake.input.regionprops_parquets]
df_combined = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
df_combined.to_parquet(snakemake.output.regionprops_aggregated_parquet, index=False)
