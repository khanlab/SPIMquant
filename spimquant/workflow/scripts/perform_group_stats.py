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

Subject filtering via ``--group-stats-where``:
An optional pandas query expression (``snakemake.params.group_stats_where``)
is applied to the merged participant dataframe *before* model fitting.  This
defines the inference cohort: subjects that do not match the expression are
excluded from both model fitting and contrast enumeration.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import logging
import os
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import dmatrix


def _setup_logging(log_path):
    """Configure logging to write to both the log file and stdout.

    Parameters
    ----------
    log_path : str or None
        Path to the Snakemake log file.  If *None* (e.g. during testing),
        only stdout is used.

    Returns
    -------
    logging.Logger
    """
    fmt = "%(asctime)s [%(name)s] [%(levelname)s] %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"
    logger = logging.getLogger("perform_group_stats")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    # Remove any handlers added by previous calls (e.g. during testing).
    logger.handlers.clear()
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.DEBUG)
    logger.addHandler(stdout_handler)
    if log_path:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w")
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)
    for handler in logger.handlers:
        handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))
    return logger


def apply_group_stats_filter(df, where_expr, log):
    """Apply a pandas query expression to filter rows from a participants DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to filter (typically the merged segstats + participant metadata).
    where_expr : str or None
        Pandas query expression string (e.g. ``"treatment in ['pbs','lecanemab']"``).
        Pass ``None`` or an empty string to skip filtering and return *df* unchanged.
    log : logging.Logger

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame (a new object; the input is not modified).

    Raises
    ------
    ValueError
        If *where_expr* is syntactically or semantically invalid, or if the
        result is empty after filtering.
    """
    if not where_expr:
        return df

    n_before = len(df)
    try:
        filtered = df.query(where_expr)
    except Exception as exc:
        raise ValueError(
            f"--group-stats-where expression is invalid: {where_expr!r}\n"
            f"Error: {exc}"
        ) from exc

    n_after = len(filtered)
    log.info(
        "[group-stats-where] Filter '%s': %d rows → %d rows (%d excluded).",
        where_expr,
        n_before,
        n_after,
        n_before - n_after,
    )

    if n_after == 0:
        raise ValueError(
            f"--group-stats-where expression '{where_expr}' excluded all rows. "
            "No subjects remain for analysis. Check the expression and "
            "the participants.tsv columns."
        )

    return filtered


def load_segstats_with_metadata(segstats_paths, participants_df, log):
    """Load all segstats files and merge with participant metadata.

    Parameters
    ----------
    segstats_paths : list
        List of paths to segstats.tsv files.
    participants_df : pd.DataFrame
        DataFrame containing participant metadata from participants.tsv.
    log : logging.Logger

    Returns
    -------
    pd.DataFrame
        Combined dataframe with segstats and participant metadata.
    """
    all_data = []

    for path in segstats_paths:
        if not os.path.exists(path):
            log.warning("Segstats file not found, skipping: %s", path)
            continue

        parts = Path(path).parts
        subject_id = next((p for p in parts if p.startswith("sub-")), None)
        if subject_id is None:
            log.warning("Could not extract subject ID from path: %s", path)
            continue

        df = pd.read_csv(path, sep="\t")
        df["participant_id"] = subject_id
        all_data.append(df)
        log.debug("Loaded segstats for %s: %d rows, %d cols", subject_id, df.shape[0], df.shape[1])

    if not all_data:
        raise ValueError("No valid segstats files found")

    combined = pd.concat(all_data, ignore_index=True)
    log.info(
        "Loaded %d subjects → combined shape: %d rows × %d cols",
        len(all_data),
        combined.shape[0],
        combined.shape[1],
    )

    merged = combined.merge(participants_df, on="participant_id", how="left")
    # Identify rows where none of the participants_df columns were matched.
    participants_cols = [c for c in participants_df.columns if c != "participant_id"]
    if participants_cols:
        unmatched_mask = merged[participants_cols].isna().all(axis=1)
        unmatched = merged.loc[unmatched_mask, "participant_id"].unique()
        if len(unmatched):
            log.warning(
                "%d participants in segstats could not be matched to participants.tsv: %s",
                len(unmatched),
                list(unmatched),
            )
    return merged


def build_prediction_row(region_data, pairwise_factor, level, strata):
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
    log,
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
    log : logging.Logger

    Returns
    -------
    dict  with keys ``{metric}_tstat``, ``{metric}_pval``, ``{metric}_cohensd``,
          ``{metric}_mean_{level_a}``, ``{metric}_mean_{level_b}``,
          ``n_{level_a}``, ``n_{level_b}``.
    """
    # Replace the placeholder with Q("column") to safely handle metric names
    # that contain special characters (e.g. '+', spaces) which patsy would
    # otherwise misinterpret as formula operators when backtick-quoted.
    actual_formula = formula_template.replace("metric", f'Q("{metric}")')

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
        log.debug(
            "  metric='%s': skipping model fit — insufficient observations "
            "(n_%s=%d, n_%s=%d, need ≥2 per group).",
            metric,
            level_a,
            n_a,
            level_b,
            n_b,
        )
        return result

    try:
        fitted = smf.ols(actual_formula, data=region_data).fit()
        log.debug(
            "  metric='%s': OLS fit OK (nobs=%d, R²=%.4f, formula=%s).",
            metric,
            fitted.nobs,
            fitted.rsquared,
            actual_formula,
        )

        # Build prediction rows for each level at the desired strata.
        pred_df_a = build_prediction_row(region_data, pairwise_factor, level_a, strata)
        pred_df_b = build_prediction_row(region_data, pairwise_factor, level_b, strata)

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

        log.debug(
            "  metric='%s': tstat=%.4f, pval=%.6f, cohensd=%.4f.",
            metric,
            tstat,
            pval,
            cohensd,
        )

        result[f"{metric}_tstat"] = tstat
        result[f"{metric}_pval"] = pval
        result[f"{metric}_cohensd"] = cohensd

    except Exception:  # noqa: BLE001
        log.warning(
            "  metric='%s': model fitting FAILED — tstat/pval will be NaN. "
            "Formula: %s  |  factor: %s  |  %s vs %s  |  strata: %s\n%s",
            metric,
            actual_formula,
            pairwise_factor,
            level_a,
            level_b,
            strata,
            traceback.format_exc(),
        )

    return result


def perform_model_based_contrast(
    data, formula, pairwise_factor, level_a, level_b, strata, metrics, log
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
    log : logging.Logger

    Returns
    -------
    pd.DataFrame
        One row per region with columns for each metric's statistics.
    """
    regions = data[["index", "name"]].drop_duplicates()
    log.info("Computing contrasts for %d regions and %d metrics.", regions.shape[0], len(metrics))
    rows = []
    n_nan_total = 0
    n_ok_total = 0

    for _, region in regions.iterrows():
        region_data = data[data["index"] == region["index"]].copy()
        log.debug("Region '%s' (index=%s): %d rows.", region["name"], region["index"], len(region_data))

        row = {"index": region["index"], "name": region["name"]}

        for metric in metrics:
            if metric not in region_data.columns:
                log.debug("  metric='%s': column not found — skipping.", metric)
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
                n_nan_total += 1
                continue

            stats = compute_contrast_for_metric(
                region_data,
                formula,
                metric,
                pairwise_factor,
                level_a,
                level_b,
                strata,
                log,
            )
            row.update(stats)
            if np.isnan(stats.get(f"{metric}_tstat", np.nan)):
                n_nan_total += 1
            else:
                n_ok_total += 1

        rows.append(row)

    log.info(
        "Finished: %d region×metric combinations succeeded, %d produced NaN stats.",
        n_ok_total,
        n_nan_total,
    )
    return pd.DataFrame(rows)


def main():
    """Main function - uses snakemake object provided by Snakemake workflow."""
    log_path = snakemake.log[0] if snakemake.log else None
    log = _setup_logging(log_path)

    log.info("=== perform_group_stats started ===")
    log.info("Output: %s", snakemake.output.stats_tsv)

    participants_df = pd.read_csv(snakemake.input.participants_tsv, sep="\t")
    log.info(
        "participants.tsv: %d rows, columns: %s",
        len(participants_df),
        list(participants_df.columns),
    )

    if "participant_id" not in participants_df.columns:
        raise ValueError("participants.tsv must contain a 'participant_id' column")

    combined_data = load_segstats_with_metadata(
        snakemake.input.segstats_tsvs, participants_df, log
    )

    formula = snakemake.params.model
    contrast_info = snakemake.params.pairwise_contrast_info
    metrics = snakemake.params.metric_columns + snakemake.params.coloc_metric_columns
    where_expr = snakemake.params.get("group_stats_where", None)

    log.info("Formula: %s", formula)
    log.info("Contrast info: %s", contrast_info)
    log.info("Metrics (%d): %s", len(metrics), metrics)
    log.info("group_stats_where: %s", where_expr)

    # Apply subject-level filter before model fitting and contrast evaluation.
    combined_data = apply_group_stats_filter(combined_data, where_expr, log)

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

    log.info(
        "Contrast: factor='%s', %s vs %s, strata=%s",
        pairwise_factor,
        level_a,
        level_b,
        strata,
    )

    # Validate required columns exist
    for col in [pairwise_factor] + list(strata.keys()):
        if col not in combined_data.columns:
            raise ValueError(
                f"Column '{col}' not found in data after merging with "
                f"participants.tsv. Available columns: {list(combined_data.columns)}"
            )

    # Log factor level distribution
    level_counts = combined_data[pairwise_factor].value_counts()
    for lvl, cnt in level_counts.items():
        log.info("  Factor '%s' level '%s': %d rows", pairwise_factor, lvl, cnt)

    # Validate that both contrast levels are still present after filtering.
    remaining_levels = set(
        str(v) for v in combined_data[pairwise_factor].dropna().unique()
    )
    for level in (level_a, level_b):
        if level not in remaining_levels:
            raise ValueError(
                f"Contrast level '{level}' for factor '{pairwise_factor}' is not "
                f"present in the filtered data. Remaining levels: {sorted(remaining_levels)}. "
                "Check --group-stats-where or ensure the requested level exists in "
                "participants.tsv."
            )

    # Log which metrics are actually present in the data
    present_metrics = [m for m in metrics if m in combined_data.columns]
    missing_metrics = [m for m in metrics if m not in combined_data.columns]
    if missing_metrics:
        log.warning(
            "%d requested metric(s) not found in segstats columns: %s",
            len(missing_metrics),
            missing_metrics,
        )
    log.info("%d / %d requested metrics are present in the data.", len(present_metrics), len(metrics))

    results = perform_model_based_contrast(
        combined_data,
        formula,
        pairwise_factor,
        level_a,
        level_b,
        strata,
        metrics,
        log,
    )

    results.to_csv(snakemake.output.stats_tsv, sep="\t", index=False)
    log.info("Wrote %d rows to %s", len(results), snakemake.output.stats_tsv)
    log.info("=== perform_group_stats complete ===")


if __name__ == "__main__":
    main()
