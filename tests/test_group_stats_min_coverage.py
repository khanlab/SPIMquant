"""Tests for --aux_min_coverage in group statistics (apply_min_coverage)."""

import importlib.util
import logging
import pathlib

import numpy as np
import pandas as pd

_SCRIPT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "spimquant"
    / "workflow"
    / "scripts"
    / "perform_group_stats.py"
)
_spec = importlib.util.spec_from_file_location("perform_group_stats", _SCRIPT)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
apply_min_coverage = _mod.apply_min_coverage

LOG = logging.getLogger("test")
METRICS = ["Abeta+density", "Abeta+fieldfrac", "coloc+density"]


def _table(coverage):
    n = len(coverage)
    return pd.DataFrame(
        {
            "index": range(n),
            "name": [f"r{i}" for i in range(n)],
            "coverage": coverage,
            "Abeta+density": np.arange(n, dtype=float) + 1,
            "Abeta+fieldfrac": np.arange(n, dtype=float) + 10,
            "participant_id": "sub-01",
        }
    )


def test_blanks_metrics_below_threshold():
    df = _table([1.0, 0.5, 0.49, np.nan])
    out = apply_min_coverage(df, METRICS, 0.5, LOG)
    keep = [True, True, False, False]
    for metric in ["Abeta+density", "Abeta+fieldfrac"]:
        assert out[metric].notna().tolist() == keep
    # region columns and coverage are left as they are, and the input unchanged
    pd.testing.assert_series_equal(out["coverage"], df["coverage"])
    assert df["Abeta+density"].notna().all()


def test_no_coverage_column_is_unchanged():
    df = _table([0.1, 0.2]).drop(columns="coverage")
    assert apply_min_coverage(df, METRICS, 0.5, LOG) is df


def test_zero_threshold_is_unchanged():
    df = _table([0.0, 0.2])
    assert apply_min_coverage(df, METRICS, 0, LOG) is df
