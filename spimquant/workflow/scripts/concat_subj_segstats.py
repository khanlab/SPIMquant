"""Concatenate subject-level segstats TSV files across all participants.

This script merges individual per-subject segstats TSV files (e.g.
mergedsegstats.tsv or colocsegstats.tsv) into a single group-level TSV,
adding a participant_id column and joining with participant metadata from
participants.tsv. The result can be exported to external statistics tools
without any further aggregation.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import os
from pathlib import Path

import pandas as pd


def extract_participant_id(path):
    """Extract participant_id from a BIDS file path.

    Parameters
    ----------
    path : str
        Path to a BIDS file.

    Returns
    -------
    str or None
        Participant ID (e.g., ``'sub-01'``) or ``None`` if not found.
    """
    for part in Path(path).parts:
        if part.startswith("sub-"):
            return part
    return None


def load_tsvs_with_metadata(tsv_paths, participants_df):
    """Load all TSV files and merge with participant metadata.

    Parameters
    ----------
    tsv_paths : list
        List of paths to subject-level TSV files.
    participants_df : pd.DataFrame
        DataFrame containing participant metadata from participants.tsv.

    Returns
    -------
    pd.DataFrame
        Combined dataframe with all subjects' data and participant metadata.
    """
    all_data = []

    for path in tsv_paths:
        if not os.path.exists(path):
            continue

        subject_id = extract_participant_id(path)
        if subject_id is None:
            continue

        df = pd.read_csv(path, sep="\t")
        df["participant_id"] = subject_id
        all_data.append(df)

    if not all_data:
        raise ValueError(
            f"No valid TSV files found. Attempted to load {len(tsv_paths)} file(s). "
            "Check that the files exist and contain valid TSV data with subject IDs "
            "in the file paths (e.g., 'sub-01')."
        )

    combined = pd.concat(all_data, ignore_index=True)
    return combined.merge(participants_df, on="participant_id", how="left")


def main():
    """Main function - uses snakemake object provided by Snakemake workflow."""
    participants_df = pd.read_csv(snakemake.input.participants_tsv, sep="\t")

    if "participant_id" not in participants_df.columns:
        raise ValueError("participants.tsv must contain a 'participant_id' column")

    combined = load_tsvs_with_metadata(snakemake.input.segstats_tsvs, participants_df)

    os.makedirs(os.path.dirname(snakemake.output.merged_tsv) or ".", exist_ok=True)
    combined.to_csv(snakemake.output.merged_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
