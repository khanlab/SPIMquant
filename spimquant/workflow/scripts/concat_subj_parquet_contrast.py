"""Concatenate subject-level parquet files filtered by contrast.

This script concatenates parquet files (regionprops or coloc) from multiple
participants, adds a participant_id column to identify each subject's data,
merges with participant metadata from participants.tsv, and filters the data
to include only rows where the specified contrast_column matches the 
contrast_value.

This is a Snakemake script that expects the `snakemake` object to be available,
which is automatically provided when executed as part of a Snakemake workflow.
"""

import os
from pathlib import Path

import pandas as pd


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


def load_parquets_with_metadata(parquet_paths, participants_df):
    """Load all parquet files and merge with participant metadata.

    Parameters
    ----------
    parquet_paths : list
        List of paths to parquet files
    participants_df : pd.DataFrame
        DataFrame containing participant metadata

    Returns
    -------
    pd.DataFrame
        Combined dataframe with all subjects' data and participant metadata
    """
    all_data = []

    for path in parquet_paths:
        if not os.path.exists(path):
            continue

        # Extract subject ID from path
        subject_id = extract_participant_id(path)

        if subject_id is None:
            continue

        # Load parquet file
        df = pd.read_parquet(path)
        df["participant_id"] = subject_id

        all_data.append(df)

    if not all_data:
        raise ValueError(
            f"No valid parquet files found. Attempted to load {len(parquet_paths)} "
            f"file(s). Check that the files exist and contain valid parquet data with "
            f"subject IDs in the file paths (e.g., 'sub-01')."
        )

    # Combine all data
    combined = pd.concat(all_data, ignore_index=True)

    # Merge with participants metadata
    merged = combined.merge(participants_df, on="participant_id", how="left")

    return merged


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

    # Load and combine all parquet files
    combined_data = load_parquets_with_metadata(
        snakemake.input.parquet_files, participants_df
    )

    # Filter data based on contrast column and value
    filtered_data = combined_data[combined_data[contrast_column] == contrast_value]

    if filtered_data.empty:
        raise ValueError(
            f"No data found for {contrast_column}={contrast_value}. "
            f"Available values in {contrast_column}: "
            f"{combined_data[contrast_column].unique().tolist()}"
        )

    # Save filtered parquet file
    filtered_data.to_parquet(snakemake.output.parquet, index=False)


if __name__ == "__main__":
    main()
