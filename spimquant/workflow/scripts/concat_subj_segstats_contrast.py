"""Concatenate subject-level segstats.tsv files filtered by contrast and compute group averages.

This script concatenates segstats.tsv files from multiple participants, adds a 
participant_id column to identify each subject's data, merges with participant 
metadata from participants.tsv, filters the data to include only rows where the 
specified contrast_column matches the contrast_value, and computes group averages
for each atlas region and metric.

This is a Snakemake script that expects the `snakemake` object to be available,
which is automatically provided when executed as part of a Snakemake workflow.
"""

import os
from pathlib import Path

import pandas as pd
import numpy as np


def extract_participant_id(path):
    """Extract participant_id from a BIDS file path.

    Parameters
    ----------
    path : str
        Path to a BIDS file

    Returns
    -------
    str or None
        Participant ID (e.g., 'sub-01') or None if not found
    """
    parts = Path(path).parts
    for part in parts:
        if part.startswith("sub-"):
            return part
    return None


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
        subject_id = extract_participant_id(path)

        if subject_id is None:
            continue

        # Load segstats
        df = pd.read_csv(path, sep="\t")
        df["participant_id"] = subject_id

        all_data.append(df)

    if not all_data:
        raise ValueError(
            f"No valid segstats files found. Attempted to load {len(segstats_paths)} "
            f"file(s). Check that the files exist and contain valid data with "
            f"subject IDs in the file paths (e.g., 'sub-01')."
        )

    # Combine all segstats
    combined = pd.concat(all_data, ignore_index=True)

    # Merge with participants metadata
    merged = combined.merge(participants_df, on="participant_id", how="left")

    return merged


def compute_group_averages(data, metric_columns):
    """Compute group averages for each atlas region and metric.

    Parameters
    ----------
    data : pd.DataFrame
        Combined dataframe with segstats data
    metric_columns : list
        List of metric column names to average

    Returns
    -------
    pd.DataFrame
        Dataframe with group averages for each region
    """
    # Group by region (index and name)
    groupby_cols = ["index", "name"]
    
    # Build aggregation dict - only include columns that exist
    agg_dict = {}
    missing_columns = []
    for col in metric_columns:
        if col in data.columns:
            agg_dict[col] = "mean"
        else:
            missing_columns.append(col)
    
    # Warn about missing columns to aid debugging
    if missing_columns:
        print(f"Warning: The following metric columns were not found in data: {missing_columns}")
        print(f"Available columns: {list(data.columns)}")
    
    if not agg_dict:
        raise ValueError(
            f"None of the specified metric columns were found in the data. "
            f"Requested: {metric_columns}. Available: {list(data.columns)}"
        )
    
    # Compute group averages
    group_avg = data.groupby(groupby_cols, as_index=False).agg(agg_dict)
    
    return group_avg


def main():
    """Main function - uses snakemake object provided by Snakemake workflow."""
    # Load participants metadata
    participants_df = pd.read_csv(snakemake.input.participants_tsv, sep="\t")

    # Validate participants.tsv has required columns
    if "participant_id" not in participants_df.columns:
        raise ValueError("participants.tsv must contain a 'participant_id' column")

    # Get contrast filtering parameters
    contrast_column = snakemake.params.contrast_column
    contrast_value = snakemake.params.contrast_value

    # Validate contrast column exists in participants.tsv
    if contrast_column not in participants_df.columns:
        raise ValueError(
            f"Contrast column '{contrast_column}' not found in participants.tsv. "
            f"Available columns: {list(participants_df.columns)}"
        )

    # Load and combine all segstats files
    combined_data = load_segstats_with_metadata(
        snakemake.input.segstats_tsvs, participants_df
    )

    # Filter data based on contrast column and value
    filtered_data = combined_data[combined_data[contrast_column] == contrast_value]

    if filtered_data.empty:
        raise ValueError(
            f"No data found for {contrast_column}={contrast_value}. "
            f"Available values in {contrast_column}: "
            f"{combined_data[contrast_column].unique().tolist()}"
        )

    # Get metric columns to average - handle None values
    metric_columns = []
    if snakemake.params.metric_columns is not None:
        metric_columns.extend(snakemake.params.metric_columns)
    if snakemake.params.coloc_metric_columns is not None:
        metric_columns.extend(snakemake.params.coloc_metric_columns)
    
    if not metric_columns:
        raise ValueError("No metric columns specified for averaging")

    # Compute group averages
    group_avg = compute_group_averages(filtered_data, metric_columns)

    # Save averaged segstats table
    group_avg.to_csv(snakemake.output.tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
