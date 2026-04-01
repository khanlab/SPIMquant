"""Per-ROI summary QC: top-region bar plots for a single subject.

Loads the merged segmentation-statistics TSV (all stains) together with
the atlas label table, then produces bar-chart visualisations of the
top brain-regions ranked by field fraction and count for each stain.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Suffixes used to identify stain-prefixed metric columns in mergedsegstats TSV.
# Columns follow the pattern "{stain}+{metric}", e.g. "Abeta+fieldfrac".
_SUFFIX_FIELDFRAC = "+fieldfrac"
_SUFFIX_COUNT = "+count"
_SUFFIX_DENSITY = "+density"


def _top_regions(df, col, n=20, ascending=False):
    """Return the top *n* rows of *df* sorted by *col*."""
    valid = df[df[col].notna() & (df[col] != 0)]
    return valid.nlargest(n, col) if not ascending else valid.nsmallest(n, col)


def _bar_plot(ax, names, values, title, xlabel, color="steelblue"):
    """Draw a horizontal bar chart on *ax*."""
    y_pos = np.arange(len(names))
    ax.barh(y_pos, values, color=color, alpha=0.8, edgecolor="white")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()  # highest value at top
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.grid(True, axis="x", alpha=0.3)


def main():
    desc = snakemake.wildcards.desc
    seg = snakemake.wildcards.seg
    template = snakemake.wildcards.template
    subject = snakemake.wildcards.subject

    stats_df = pd.read_csv(snakemake.input.segstats, sep="\t")
    label_df = pd.read_csv(snakemake.input.label_tsv, sep="\t")

    # Drop background (atlas label 0) — those voxels are outside the brain
    if "index" in stats_df.columns:
        stats_df = stats_df[stats_df["index"] != 0].copy()

    # Merge region names (label_df has 'index' and 'name' columns)
    if "name" not in stats_df.columns and "index" in stats_df.columns:
        stats_df = stats_df.merge(label_df[["index", "name"]], on="index", how="left")

    region_name_col = "name" if "name" in stats_df.columns else "index"

    # Identify stain-prefixed metric columns (pattern: "{stain}+{metric}")
    ff_cols = [c for c in stats_df.columns if c.endswith(_SUFFIX_FIELDFRAC)]
    count_cols = [c for c in stats_df.columns if c.endswith(_SUFFIX_COUNT)]
    density_cols = [c for c in stats_df.columns if c.endswith(_SUFFIX_DENSITY)]

    # Determine number of rows: 1 row per metric type (ff, count, density)
    # with one subplot per stain within each row
    n_ff = len(ff_cols)
    n_count = len(count_cols)
    n_density = len(density_cols)
    n_rows = (1 if n_ff else 0) + (1 if n_count else 0) + (1 if n_density else 0)

    if n_rows == 0:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(
            0.5,
            0.5,
            "No stain-prefixed metric columns found in segstats TSV",
            ha="center",
            va="center",
            fontsize=12,
            color="gray",
            transform=ax.transAxes,
        )
        ax.axis("off")
        plt.savefig(snakemake.output.png, dpi=150, bbox_inches="tight")
        plt.close()
        return

    n_top = 20
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    row_specs = []
    if n_ff:
        row_specs.append(("Field Fraction (%)", ff_cols, "steelblue"))
    if n_count:
        row_specs.append(("Count (objects)", count_cols, "darkorange"))
    if n_density:
        row_specs.append(("Density (objects/vol)", density_cols, "forestgreen"))

    max_stains = max(len(cols) for _, cols, _ in row_specs)
    fig_width = max(10, max_stains * 6)
    fig_height = n_rows * 6

    fig, axes = plt.subplots(
        n_rows, max_stains, figsize=(fig_width, fig_height), squeeze=False
    )
    fig.suptitle(
        f"Per-ROI Summary QC\n"
        f"Subject: {subject}  |  Atlas: {seg}  |  Template: {template}  |  "
        f"Method: {desc}",
        fontsize=12,
        fontweight="bold",
    )

    for row_idx, (metric_label, metric_cols, base_color) in enumerate(row_specs):
        for col_idx in range(max_stains):
            ax = axes[row_idx, col_idx]
            if col_idx >= len(metric_cols):
                ax.axis("off")
                continue
            col = metric_cols[col_idx]
            stain_label = col.split("+")[0]
            top = _top_regions(stats_df, col, n=n_top)
            names = top[region_name_col].astype(str).tolist()
            values = top[col].values
            color = colors[col_idx % len(colors)]
            _bar_plot(
                ax,
                names,
                values,
                title=f"Top {n_top} regions  —  {stain_label}",
                xlabel=metric_label,
                color=color,
            )

    plt.tight_layout()
    plt.savefig(snakemake.output.png, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved per-ROI summary QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
