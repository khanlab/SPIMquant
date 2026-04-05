"""Perform formula-based group statistical analysis on segmentation statistics.

This script reads segstats.tsv files from multiple participants, fits a single
global OLS model per region/metric using the user-supplied patsy/statsmodels
formula, and computes pairwise contrast statistics using the model's covariance
matrix.

The contrast is specified via ``snakemake.params.pairwise_contrast_info``, a
dict with keys:
- ``factor``: the column in participants.tsv to compare levels of
- ``levelA`` / ``levelB``: the two levels to contrast (levelA − levelB)
- ``strata``: a dict of {column: value} pairs for stratified analyses;
  if non-empty the model is fit with *all* data but marginal means are
  evaluated at the given stratum values.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrix


def load_segstats_with_metadata(segstats_paths, participants_df):
    """Load all segstats files and merge with participant metadata.

    Parameters
    ----------
    segstats_paths : list
        List of paths to segstats.tsv files.
    participants_df : pd.DataFrame
        DataFrame containing participant metadata from participants.tsv.

    Returns
    -------
    pd.DataFrame
        Combined dataframe with segstats and participant metadata.
    """
    all_data = []

    for path in segstats_paths:
        if not os.path.exists(path):
            continue

        parts = Path(path).parts
        subject_id = next((p for p in parts if p.startswith("sub-")), None)
        if subject_id is None:
            continue

        df = pd.read_csv(path, sep="\t")
        df["participant_id"] = subject_id
        all_data.append(df)

    if not all_data:
        raise ValueError("No valid segstats files found")

    combined = pd.concat(all_data, ignore_index=True)
    return combined.merge(participants_df, on="participant_id", how="left")


def _build_prediction_row(region_data, pairwise_factor, level, strata):
    """Build a one-row DataFrame for marginal mean prediction.

    Continuous variables are held at their mean; categorical variables are
    held at their mode.  The pairwise factor and any strata variables are set
    to the supplied values.

    Parameters
    ----------
    region_data : pd.DataFrame
    pairwise_factor : str
    level : str
    strata : dict

    Returns
    -------
    pd.DataFrame  (single row)
    """
    row = {}
    for col in region_data.columns:
        if col == "participant_id":
            continue
        if pd.api.types.is_numeric_dtype(region_data[col]):
            row[col] = [region_data[col].mean()]
        else:
            mode_vals = region_data[col].mode()
            row[col] = [mode_vals.iloc[0] if len(mode_vals) > 0 else None]

    row[pairwise_factor] = [level]
    for factor, value in strata.items():
        row[factor] = [value]

    return pd.DataFrame(row)


def compute_contrast_for_metric(
    region_data,
    formula_template,
    metric,
    pairwise_factor,
    level_a,
    level_b,
    strata,
):
    """Fit a global OLS model and compute one pairwise contrast for *metric*.

    The model is fit on all rows of *region_data*.  Marginal means are
    computed by predicting at reference covariate values with the pairwise
    factor set to each level in turn (and strata variables fixed at the
    requested values).

    Parameters
    ----------
    region_data : pd.DataFrame
    formula_template : str
        Formula with the literal string ``metric`` as the response variable
        placeholder, e.g. ``"metric ~ C(treatment) + age"``.
    metric : str
        Actual metric column name; replaces ``metric`` in the formula.
    pairwise_factor : str
    level_a, level_b : str
    strata : dict

    Returns
    -------
    dict  with keys ``{metric}_tstat``, ``{metric}_pval``, ``{metric}_cohensd``,
          ``{metric}_mean_{level_a}``, ``{metric}_mean_{level_b}``,
          ``n_{level_a}``, ``n_{level_b}``.
    """
    # Replace the placeholder with the actual (possibly backtick-quoted) column.
    actual_formula = formula_template.replace("metric", f"`{metric}`")

    # Filter to rows relevant for the raw-mean / Cohen's d calculation.
    grp_a = region_data[region_data[pairwise_factor] == level_a]
    grp_b = region_data[region_data[pairwise_factor] == level_b]
    if strata:
        for f, v in strata.items():
            grp_a = grp_a[grp_a[f] == v]
            grp_b = grp_b[grp_b[f] == v]

    n_a = int(grp_a[metric].dropna().shape[0])
    n_b = int(grp_b[metric].dropna().shape[0])
    mean_a = float(grp_a[metric].mean()) if n_a > 0 else np.nan
    mean_b = float(grp_b[metric].mean()) if n_b > 0 else np.nan
    std_a = float(grp_a[metric].std()) if n_a > 1 else np.nan
    std_b = float(grp_b[metric].std()) if n_b > 1 else np.nan

    result = {
        f"n_{level_a}": n_a,
        f"n_{level_b}": n_b,
        f"{metric}_mean_{level_a}": mean_a,
        f"{metric}_mean_{level_b}": mean_b,
        f"{metric}_tstat": np.nan,
        f"{metric}_pval": np.nan,
        f"{metric}_cohensd": np.nan,
    }

    if n_a < 2 or n_b < 2:
        return result

    try:
        fitted = smf.ols(actual_formula, data=region_data).fit()

        # Build prediction rows for each level at the desired strata.
        pred_df_a = _build_prediction_row(region_data, pairwise_factor, level_a, strata)
        pred_df_b = _build_prediction_row(region_data, pairwise_factor, level_b, strata)

        # Use patsy with the model's design_info for consistent dummy encoding.
        design_info = fitted.model.data.design_info
        dm_a = np.asarray(dmatrix(design_info, pred_df_a, return_type="matrix"))
        dm_b = np.asarray(dmatrix(design_info, pred_df_b, return_type="matrix"))
        contrast_vec = (dm_a - dm_b)[0]

        ct = fitted.t_test(contrast_vec)
        tstat = float(np.asarray(ct.tvalue).item())
        pval = float(np.asarray(ct.pvalue).item())

        # Cohen's d from pooled standard deviation.
        denom = n_a + n_b - 2
        if denom > 0 and not (np.isnan(std_a) or np.isnan(std_b)):
            pooled_var = ((n_a - 1) * std_a**2 + (n_b - 1) * std_b**2) / denom
            pooled_std = np.sqrt(pooled_var) if pooled_var >= 0 else np.nan
            cohensd = (mean_a - mean_b) / pooled_std if pooled_std > 0 else np.nan
        else:
            cohensd = np.nan

        result[f"{metric}_tstat"] = tstat
        result[f"{metric}_pval"] = pval
        result[f"{metric}_cohensd"] = cohensd

    except Exception:  # noqa: BLE001
        pass  # leave NaN placeholders

    return result


def perform_model_based_contrast(
    data, formula, pairwise_factor, level_a, level_b, strata, metrics
):
    """Compute pairwise contrast statistics for every region and metric.

    Parameters
    ----------
    data : pd.DataFrame
        Combined dataframe with segstats and participant metadata.
    formula : str
        Model formula (patsy/statsmodels), with ``metric`` as placeholder.
    pairwise_factor : str
    level_a, level_b : str
    strata : dict
    metrics : list[str]

    Returns
    -------
    pd.DataFrame
        One row per region with columns for each metric's statistics.
    """
    regions = data[["index", "name"]].drop_duplicates()
    rows = []

    for _, region in regions.iterrows():
        region_data = data[data["index"] == region["index"]].copy()

        row = {"index": region["index"], "name": region["name"]}

        for metric in metrics:
            if metric not in region_data.columns:
                row.update(
                    {
                        f"n_{level_a}": np.nan,
                        f"n_{level_b}": np.nan,
                        f"{metric}_mean_{level_a}": np.nan,
                        f"{metric}_mean_{level_b}": np.nan,
                        f"{metric}_tstat": np.nan,
                        f"{metric}_pval": np.nan,
                        f"{metric}_cohensd": np.nan,
                    }
                )
                continue

            stats = compute_contrast_for_metric(
                region_data,
                formula,
                metric,
                pairwise_factor,
                level_a,
                level_b,
                strata,
            )
            row.update(stats)

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    """Main function - uses snakemake object provided by Snakemake workflow."""
    participants_df = pd.read_csv(snakemake.input.participants_tsv, sep="\t")

    if "participant_id" not in participants_df.columns:
        raise ValueError("participants.tsv must contain a 'participant_id' column")

    combined_data = load_segstats_with_metadata(
        snakemake.input.segstats_tsvs, participants_df
    )

    formula = snakemake.params.model
    contrast_info = snakemake.params.pairwise_contrast_info
    metrics = snakemake.params.metric_columns + snakemake.params.coloc_metric_columns

    if not formula:
        raise ValueError("No model formula supplied (--model is required).")
    if not contrast_info:
        raise ValueError(
            "No pairwise contrast information found for this wildcard. "
            "Check that --pairwise matches a column in participants.tsv."
        )

    pairwise_factor = contrast_info["factor"]
    level_a = contrast_info["levelA"]
    level_b = contrast_info["levelB"]
    strata = contrast_info.get("strata", {})

    # Validate required columns exist
    for col in [pairwise_factor] + list(strata.keys()):
        if col not in combined_data.columns:
            raise ValueError(
                f"Column '{col}' not found in data after merging with "
                f"participants.tsv. Available columns: {list(combined_data.columns)}"
            )

    results = perform_model_based_contrast(
        combined_data,
        formula,
        pairwise_factor,
        level_a,
        level_b,
        strata,
        metrics,
    )

    results.to_csv(snakemake.output.stats_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
