"""Tests for --group-stats-where filtering in perform_group_stats.py.

Covers:
- apply_group_stats_filter passes through data unchanged when no expression given
- apply_group_stats_filter reduces rows matching a valid expression
- apply_group_stats_filter raises ValueError on an invalid expression
- apply_group_stats_filter raises ValueError when all rows are excluded
- Planning-time enumeration uses the filtered participants dataframe so that
  contrast labels only cover the requested levels
- Planning-time enumeration raises ValueError when fewer than 2 levels remain
"""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pytest


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


REPO_ROOT = _find_repo_root(Path(__file__).parent)
SCRIPT_PATH = REPO_ROOT / "spimquant/workflow/scripts/perform_group_stats.py"

_spec = spec_from_file_location("perform_group_stats", SCRIPT_PATH)
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)

apply_group_stats_filter = _mod.apply_group_stats_filter


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def sample_df():
    """A small participants-style dataframe with a treatment column."""
    return pd.DataFrame(
        {
            "participant_id": [
                "sub-01",
                "sub-02",
                "sub-03",
                "sub-04",
                "sub-05",
                "sub-06",
            ],
            "treatment": ["pbs", "pbs", "lecanemab", "lecanemab", "control", "n/a"],
            "age": [12, 13, 11, 14, 12, 15],
        }
    )


# ---------------------------------------------------------------------------
# apply_group_stats_filter unit tests
# ---------------------------------------------------------------------------


def test_filter_none_returns_unchanged(sample_df):
    """Passing None as expression returns the dataframe unchanged."""
    result = apply_group_stats_filter(sample_df, None)
    pd.testing.assert_frame_equal(result, sample_df)


def test_filter_empty_string_returns_unchanged(sample_df):
    """Passing an empty string returns the dataframe unchanged."""
    result = apply_group_stats_filter(sample_df, "")
    pd.testing.assert_frame_equal(result, sample_df)


def test_filter_reduces_rows(sample_df):
    """A valid filter expression keeps only matching rows."""
    result = apply_group_stats_filter(
        sample_df, "treatment in ['pbs', 'lecanemab']"
    )
    assert len(result) == 4
    assert set(result["treatment"].unique()) == {"pbs", "lecanemab"}


def test_filter_single_condition(sample_df):
    """A simple equality filter keeps only matching rows."""
    result = apply_group_stats_filter(sample_df, "age >= 13")
    assert all(result["age"] >= 13)


def test_filter_invalid_expression_raises(sample_df):
    """An invalid/unparseable expression raises a ValueError with a clear message."""
    with pytest.raises(ValueError, match="--group-stats-where expression is invalid"):
        apply_group_stats_filter(sample_df, "treatment @@@ 'pbs'")


def test_filter_nonexistent_column_raises(sample_df):
    """Referencing a column that doesn't exist raises a ValueError."""
    with pytest.raises(ValueError, match="--group-stats-where expression is invalid"):
        apply_group_stats_filter(sample_df, "nonexistent_col == 'pbs'")


def test_filter_all_excluded_raises(sample_df):
    """Filtering to zero rows raises a ValueError with a clear message."""
    with pytest.raises(ValueError, match="excluded all rows"):
        apply_group_stats_filter(sample_df, "treatment == 'does_not_exist'")


def test_filter_preserves_row_order(sample_df):
    """Filtered result preserves the relative ordering of rows."""
    result = apply_group_stats_filter(sample_df, "age < 14")
    # Rows 0, 1, 2, 4 have age < 14
    expected_ids = ["sub-01", "sub-02", "sub-03", "sub-05"]
    assert list(result["participant_id"]) == expected_ids


# ---------------------------------------------------------------------------
# Planning-time contrast enumeration tests (Snakefile logic extracted)
# ---------------------------------------------------------------------------
# The Snakefile logic is not importable as a module, so we reproduce the
# essential filtering + enumeration logic here using the same utility function
# to verify the end-to-end behaviour.


from itertools import combinations as _combinations


def _enumerate_pairwise_contrasts(participants_df, pairwise_factors, where_expr=None):
    """Simplified version of the planning-time contrast enumeration in Snakefile."""
    df = apply_group_stats_filter(participants_df, where_expr)

    labels = []
    info = {}
    for factor in pairwise_factors:
        if factor not in df.columns:
            raise ValueError(f"Pairwise factor '{factor}' not found.")
        levels = sorted(str(v) for v in df[factor].dropna().unique())
        if len(levels) < 2:
            raise ValueError(
                f"Pairwise factor '{factor}' has fewer than 2 levels after "
                f"applying --group-stats-where filter. Remaining: {levels}"
            )
        for lA, lB in _combinations(levels, 2):
            label = f"{factor}+{lA}vs{lB}"
            labels.append(label)
            info[label] = {"factor": factor, "levelA": lA, "levelB": lB, "strata": {}}
    return labels, info


def test_planning_no_filter_all_contrasts(sample_df):
    """Without a filter, all unique level pairs are enumerated."""
    labels, _ = _enumerate_pairwise_contrasts(sample_df, ["treatment"])
    # treatment has 4 levels → C(4,2) = 6 contrasts
    assert len(labels) == 6


def test_planning_with_filter_fewer_contrasts(sample_df):
    """With a filter keeping only pbs/lecanemab, only 1 contrast is enumerated."""
    labels, info = _enumerate_pairwise_contrasts(
        sample_df, ["treatment"], "treatment in ['pbs', 'lecanemab']"
    )
    assert len(labels) == 1
    # Both contrast levels must be pbs and lecanemab (in sorted order)
    label = labels[0]
    assert label.startswith("treatment+")
    assert "pbs" in label and "lecanemab" in label


def test_planning_filter_fewer_than_two_levels_raises(sample_df):
    """Filter leaving < 2 levels for a pairwise factor raises ValueError."""
    with pytest.raises(ValueError, match="fewer than 2 levels"):
        _enumerate_pairwise_contrasts(
            sample_df, ["treatment"], "treatment == 'pbs'"
        )


def test_planning_invalid_filter_raises(sample_df):
    """An invalid filter expression during planning raises ValueError."""
    with pytest.raises(ValueError, match="--group-stats-where expression is invalid"):
        _enumerate_pairwise_contrasts(
            sample_df, ["treatment"], "@@invalid@@"
        )
