import pandas as pd

# read each TSV into a list of dataframes
dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input]

# merge all dataframes on 'index'
from functools import reduce

merged = reduce(lambda left, right: pd.merge(left, right, on=("index", "name"), how="outer"), dfs)

# fill all columns except volume with 0
cols_to_fill = [c for c in merged.columns if c != "volume"]
merged[cols_to_fill] = merged[cols_to_fill].fillna(0)

# compute density only if columns exist
if "count" in merged.columns and "volume" in merged.columns:
    merged["density"] = merged["count"] / merged["volume"]


merged.to_csv(snakemake.output.tsv, sep='\t',index=False)
