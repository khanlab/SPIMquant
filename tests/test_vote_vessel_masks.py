from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import dask.array as da
import numpy as np


def _find_repo_root(start: Path) -> Path:
    """Find repository root by locating pyproject.toml."""
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


REPO_ROOT = _find_repo_root(Path(__file__).parent)
SCRIPT_PATH = REPO_ROOT / "spimquant/workflow/scripts/vote_vessel_masks.py"


# This script is executed by Snakemake as a standalone file, not as an importable
# package module; we use importlib to unit-test helper functions directly.
spec = spec_from_file_location("vote_vessel_masks", SCRIPT_PATH)
vote_mod = module_from_spec(spec)
spec.loader.exec_module(vote_mod)


def _masks(*rows):
    """Build (1, 1, 1, N) 0/100 masks from strings like '1100'."""
    return [
        (np.array([int(c) for c in r], dtype=np.uint8) * 100).reshape(1, 1, 1, -1)
        for r in rows
    ]


def _vote(arrays):
    return np.asarray(vote_mod.majority_vote(arrays).compute()).ravel()


def test_three_masks_need_two_votes():
    masks = _masks("11100000", "10011000", "01010100")
    # votes:          2 2 1 2 1 1 0 0
    np.testing.assert_array_equal(_vote(masks), [100, 100, 0, 100, 0, 0, 0, 0])


def test_two_masks_give_intersection():
    masks = _masks("1100", "1010")
    np.testing.assert_array_equal(_vote(masks), [100, 0, 0, 0])


def test_any_positive_value_counts_as_foreground():
    a = np.array([1, 0, 255], dtype=np.uint8)
    b = np.array([100, 0, 1], dtype=np.uint8)
    c = np.zeros(3, dtype=np.uint8)
    np.testing.assert_array_equal(_vote([a, b, c]), [100, 0, 100])


def test_output_is_uint8_zero_or_hundred():
    rng = np.random.default_rng(0)
    masks = [(rng.random((1, 8, 8, 8)) > 0.5).astype(np.uint8) * 100 for _ in range(3)]
    out = vote_mod.majority_vote(masks).compute()
    assert out.dtype == np.uint8
    assert set(np.unique(out)) <= {0, 100}
    assert out.shape == masks[0].shape


def test_dask_and_numpy_inputs_agree():
    rng = np.random.default_rng(1)
    masks = [(rng.random((1, 16, 16, 16)) > 0.4).astype(np.uint8) * 100 for _ in range(3)]
    dmasks = [da.from_array(m, chunks=(1, 8, 8, 8)) for m in masks]
    np.testing.assert_array_equal(_vote(masks), _vote(dmasks))
