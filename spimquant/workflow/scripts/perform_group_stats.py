"""Perform group-based statistical analysis on segmentation statistics.

This script reads segstats.tsv files from multiple participants and performs
statistical tests (e.g., t-tests, ANOVA) based on contrasts defined in the
participants.tsv file.

This is a Snakemake script that expects the `snakemake` object to be available,
which is automatically provided when executed as part of a Snakemake workflow.
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path


def load_segstats_with_metadata(segstats_paths, participants_df):
    """Load all segstats files and merge with participant metadata.

    Parameters
    ----------
    segstats_paths : list
        List of paths to segstats.tsv files
    participants_df : pd.DataFrame
        DataFrame containing participant metadata

    Returns
    -------
    pd.DataFrame
        Combined dataframe with segstats and participant metadata
    """
    all_data = []

    for path in segstats_paths:
        if not os.path.exists(path):
            continue

        # Extract subject ID from path
        # Path format: ...sub-{subject_id}/...
        parts = Path(path).parts
        subject_id = None
        for part in parts:
            if part.startswith("sub-"):
                subject_id = part
                break

        if subject_id is None:
            continue

        # Load segstats
        df = pd.read_csv(path, sep="\t")
        df["participant_id"] = subject_id

        all_data.append(df)

    if not all_data:
        raise ValueError("No valid segstats files found")

    # Combine all segstats
    combined = pd.concat(all_data, ignore_index=True)

    # Merge with participants metadata
    merged = combined.merge(participants_df, on="participant_id", how="left")

    return merged


def perform_two_group_test(data, group_column, group1_value, group2_value, metrics):
    """Perform two-sample t-tests for each region and metric.

    Parameters
    ----------
    data : pd.DataFrame
        Combined dataframe with segstats and metadata
    group_column : str
        Column name for grouping (e.g., 'treatment')
    group1_value : str
        Value for group 1 (e.g., 'control')
    group2_value : str
        Value for group 2 (e.g., 'drug')
    metrics : list
        List of metric columns to test (e.g., ['fieldfrac', 'density'])

    Returns
    -------
    pd.DataFrame
        Results dataframe with statistics for each region
    """
    results = []

    # Get unique regions
    regions = data[["index", "name"]].drop_duplicates()

    for _, region in regions.iterrows():
        region_idx = region["index"]
        region_name = region["name"]

        # Filter data for this region
        region_data = data[data["index"] == region_idx]

        result = {
            "index": region_idx,
            "name": region_name,
        }

        # Get group data
        group1_data = region_data[region_data[group_column] == group1_value]
        group2_data = region_data[region_data[group_column] == group2_value]

        result[f"n_{group1_value}"] = len(group1_data)
        result[f"n_{group2_value}"] = len(group2_data)

        # Perform tests for each metric
        for metric in metrics:
            if metric not in region_data.columns:
                continue

            g1_values = group1_data[metric].dropna()
            g2_values = group2_data[metric].dropna()

            if len(g1_values) < 2 or len(g2_values) < 2:
                # Not enough data for testing
                result[f"{metric}_mean_{group1_value}"] = (
                    g1_values.mean() if len(g1_values) > 0 else np.nan
                )
                result[f"{metric}_mean_{group2_value}"] = (
                    g2_values.mean() if len(g2_values) > 0 else np.nan
                )
                result[f"{metric}_tstat"] = np.nan
                result[f"{metric}_pval"] = np.nan
                continue

            # Calculate means
            result[f"{metric}_mean_{group1_value}"] = g1_values.mean()
            result[f"{metric}_mean_{group2_value}"] = g2_values.mean()
            result[f"{metric}_std_{group1_value}"] = g1_values.std()
            result[f"{metric}_std_{group2_value}"] = g2_values.std()

            # Perform two-sample t-test
            tstat, pval = stats.ttest_ind(g1_values, g2_values)
            result[f"{metric}_tstat"] = tstat
            result[f"{metric}_pval"] = pval

            # Calculate effect size (Cohen's d)
            pooled_std = np.sqrt(
                (
                    (len(g1_values) - 1) * g1_values.std() ** 2
                    + (len(g2_values) - 1) * g2_values.std() ** 2
                )
                / (len(g1_values) + len(g2_values) - 2)
            )
            if pooled_std > 0:
                cohens_d = (g1_values.mean() - g2_values.mean()) / pooled_std
                result[f"{metric}_cohensd"] = cohens_d
            else:
                result[f"{metric}_cohensd"] = np.nan

        results.append(result)

    return pd.DataFrame(results)


def main():
    """Main function - uses snakemake object provided by Snakemake workflow."""
    participants_df = pd.read_csv(snakemake.input.participants_tsv, sep="\t")

    # Validate participants.tsv has required columns
    if "participant_id" not in participants_df.columns:
        raise ValueError("participants.tsv must contain a 'participant_id' column")

    # Load and combine all segstats files
    combined_data = load_segstats_with_metadata(
        snakemake.input.segstats_tsvs, participants_df
    )

    # Get contrast information
    contrast_column = snakemake.params.contrast_column
    contrast_values = snakemake.params.contrast_values
    metrics = snakemake.params.metric_columns + snakemake.params.coloc_metric_columns

    # Validate contrast column exists if specified
    if contrast_column is not None and contrast_column not in participants_df.columns:
        raise ValueError(
            f"Contrast column '{contrast_column}' not found in participants.tsv. "
            f"Available columns: {list(participants_df.columns)}"
        )

    if contrast_column is None or contrast_values is None:
        # No contrasts specified - just aggregate data
        # Group by region and compute summary statistics

        agg_dict = {"participant_id": "count"}
        for metric in metrics:
            agg_dict[metric] = ["mean", "std", "min", "max"]

        results = combined_data.groupby(["index", "name"]).agg(agg_dict).reset_index()
        results.columns = ["_".join(col).strip("_") for col in results.columns.values]
        results.rename(columns={"participant_id_count": "n_subjects"}, inplace=True)

    elif len(contrast_values) == 2:
        # Two-group comparison

        results = perform_two_group_test(
            combined_data,
            contrast_column,
            contrast_values[0],
            contrast_values[1],
            metrics,
        )
    else:
        raise ValueError(
            "Currently only two-group contrasts are supported. "
            f"Got {len(contrast_values)} groups: {contrast_values}"
        )

    # Save results
    results.to_csv(snakemake.output.stats_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
