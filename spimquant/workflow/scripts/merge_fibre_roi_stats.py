"""Merge per-stain fibre ROI statistics into one segstats-style table.

Columns are ``index``, ``name`` and ``{stain}+{metric}`` (the naming the
group statistics, heatmap and featuremap rules use for segmentation
statistics), so the result can be written as ``desc-frangi_mergedsegstats``
and go through the same group-level rules.
"""

from functools import reduce

import pandas as pd


def merge_tables(dfs, stains):
    renamed = []
    for df, stain in zip(dfs, stains):
        cols = {c: f"{stain}+{c}" for c in df.columns if c not in ("index", "name")}
        renamed.append(df.rename(columns=cols))
    return reduce(
        lambda a, b: pd.merge(a, b, on=["index", "name"], how="outer"), renamed
    )


def main():
    dfs = [pd.read_csv(f, sep="\t") for f in snakemake.input.tsvs]
    merge_tables(dfs, snakemake.params.stains).to_csv(
        snakemake.output.tsv, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
