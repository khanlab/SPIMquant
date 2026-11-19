"""Combine per-scan landmark TRE tables; median/IQR per path x method x stain x inversion
(and per landmark). Rows with points outside a transform's field are kept in the
combined table but left out of the summaries (their TRE is not meaningful)."""
import pandas as pd

sm = snakemake  # noqa: F821
df = pd.concat([pd.read_csv(f, sep="\t") for f in sm.input.tsvs], ignore_index=True)
df.to_csv(sm.output.all_tsv, sep="\t", index=False)

ok = df[df["inside_field"].astype(bool)]
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
summ["n_outside_field"] = (
    df[~df["inside_field"].astype(bool)].groupby(keys).size().reindex(pd.MultiIndex.from_frame(summ[keys])).fillna(0).astype(int).values
)
summ.sort_values(["path", "inversion", "median_um"]).round(2).to_csv(sm.output.summary, sep="\t", index=False)

per = ok.groupby(keys + ["landmark"]).apply(stats, include_groups=False).reset_index()
per.round(2).to_csv(sm.output.per_landmark, sep="\t", index=False)
