"""Jacobian determinants in the registration QC reports (regqc_lib.py)."""

import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

pytest.importorskip("PIL")
sys.path.insert(
    0, str(Path(__file__).parents[1] / "spimquant" / "workflow" / "scripts")
)
from regqc_lib import _display_window, jacobian_determinant  # noqa: E402


def _field(affine, scale_lps, shape=(20, 24, 28)):
    """ITK-style field for x -> scale * x, with displacements in LPS mm."""
    ijk = np.stack(np.indices(shape), axis=-1).astype(np.float32)
    ras = ijk @ affine[:3, :3].T + affine[:3, 3]
    lps = ras * np.array([-1.0, -1.0, 1.0])
    disp = lps * (np.array(scale_lps) - 1)
    return nib.Nifti1Image(disp[:, :, :, None, :].astype(np.float32), affine)


@pytest.mark.parametrize(
    "affine",
    [
        np.diag([0.1, 0.1, 0.1, 1.0]),
        np.diag([-0.1, 0.2, 0.05, 1.0]),
    ],
)
def test_uniform_scaling(affine):
    det = jacobian_determinant(_field(affine, (1.1, 0.9, 1.2)), stride=2)
    assert det == pytest.approx(1.1 * 0.9 * 1.2, rel=1e-3)


def test_reflection_is_folded():
    det = jacobian_determinant(_field(np.diag([0.1, 0.1, 0.1, 1.0]), (-1, 1, 1)))
    assert (det < 0).all()


def test_canonical_orientation():
    """The determinant is reoriented like the images shown in the report."""
    affine = np.diag([-0.1, 0.1, 0.1, 1.0])
    det = jacobian_determinant(_field(affine, (1, 1, 1)), stride=2)
    assert det.shape == (10, 12, 14)
    assert det == pytest.approx(1.0, rel=1e-4)


def test_display_window_keeps_skewed_stain_visible():
    # most of the brain is dim, with a few very bright voxels (e.g. YoPro)
    rng = np.random.default_rng(0)
    data = rng.normal(0.15, 0.005, (40, 40, 40)).astype(np.float32)
    data[:5] = 0.7
    mask = np.ones(data.shape, dtype=bool)
    low, high, gamma = _display_window(data, mask)
    assert high > 0.5
    median = ((np.median(data) - low) / (high - low)) ** gamma
    assert 0.3 < median < 0.5


def test_display_window_handles_empty_image():
    assert _display_window(np.zeros((4, 4, 4)), None) == (0.0, 1.0, 1.0)
