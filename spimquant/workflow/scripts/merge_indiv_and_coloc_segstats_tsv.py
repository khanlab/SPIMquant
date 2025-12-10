
"""Merge individual stain statistics with colocalization statistics into a single TSV.

This script combines per-stain segmentation statistics (e.g., fieldfrac, density) 
with colocalization statistics (e.g., overlap_ratio, distance) into a unified TSV
file where columns are prefixed by stain name or 'coloc'.

This is a Snakemake script that expects the `snakemake` object to be available.
"""

import pandas as pd

indiv_files = snakemake.input.indiv_tsvs
coloc_file = snakemake.input.coloc_tsv
output_file = snakemake.output.merged_tsv
stains = snakemake.params.stains  # list aligned to indiv_files


# Helper: load TSV + prefix non-key columns
def load_and_prefix(path, prefix):
    df = pd.read_csv(path, sep="\t")
    # prefix all columns except index and 'name'
    rename_cols = {
        col: f"{prefix}+{col}" for col in df.columns if col not in ["index", "name"]
    }
    return df.rename(columns=rename_cols)


# Load and merge all indiv tsvs
merged = None

for stain, file in zip(stains, indiv_files):
    df = load_and_prefix(file, stain)

    # merge incrementally on ['index','name']
    if merged is None:
        merged = df
    else:
        merged = merged.merge(df, on=["index", "name"], how="outer")


# Load coloc + merge with prefix
coloc_df = load_and_prefix(coloc_file, "coloc")
merged = merged.merge(coloc_df, on=["index", "name"], how="outer")


# Save output
merged.to_csv(output_file, sep="\t", index=False)
