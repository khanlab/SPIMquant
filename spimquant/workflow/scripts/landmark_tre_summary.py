"""Combine per-scan landmark TRE tables; median/IQR per path x method x stain x inversion
(and per landmark). Rows with points outside a transform's field, and fixed-point
rows whose inverse did not converge (residual above landmark_files.INV_TOL_UM), are
kept in the combined table but left out of the summaries (their TRE is not
meaningful); the summary counts both."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import landmark_files as LF  # noqa: E402

sm = snakemake  # noqa: F821
tables = [pd.read_csv(f, sep="\t") for f in sm.input.tsvs]
df = pd.concat([t for t in tables if not t.empty] or tables, ignore_index=True)
df.to_csv(sm.output.all_tsv, sep="\t", index=False)

inside = df["inside_field"].astype(bool)
conv = LF.inv_converged(df)
ok = df[inside & conv]
keys = ["path", "method", "stain", "inversion"]


def stats(g):
    t = g["tre_um"]
    return pd.Series(
        dict(
            n=len(t),
            n_subjects=g["subject"].nunique(),
            median_um=t.median(),
            q25_um=t.quantile(0.25),
            q75_um=t.quantile(0.75),
            mean_um=t.mean(),
            sd_um=t.std(),
            max_um=t.max(),
        )
    )


summ = ok.groupby(keys).apply(stats, include_groups=False).reset_index()


def count(mask):
    idx = pd.MultiIndex.from_frame(summ[keys])
    return df[mask].groupby(keys).size().reindex(idx).fillna(0).astype(int).values


summ["n_outside_field"] = count(~inside)
summ["n_not_converged"] = count(inside & ~conv)
summ.sort_values(["path", "inversion", "median_um"]).round(2).to_csv(sm.output.summary, sep="\t", index=False)

per = ok.groupby(keys + ["landmark"]).apply(stats, include_groups=False).reset_index()
per.round(2).to_csv(sm.output.per_landmark, sep="\t", index=False)
