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
from regqc_lib import (  # noqa: E402
    _display_window,
    build_viewer_data,
    jacobian_determinant,
    render_report,
)


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


def _write(path, data):
    nib.save(nib.Nifti1Image(data.astype(np.float32), np.eye(4)), path)
    return str(path)


def _images(tmp_path):
    rng = np.random.default_rng(0)
    return [
        _write(tmp_path / f"{name}.nii.gz", rng.random((24, 26, 28)) + 0.1)
        for name in ["fixed_a", "fixed_b", "moved_a", "moved_b"]
    ]


def test_report_renders_a_moving_image_once_when_stains_share_it(tmp_path):
    fixed_a, fixed_b, moved, _ = _images(tmp_path)
    viewer = build_viewer_data(
        fixed_a,
        [{"name": "m", "warped": moved}],
        [
            {"name": "A", "fixed": fixed_a, "warped": [moved]},
            {"name": "B", "fixed": fixed_b, "warped": [moved]},
        ],
    )
    assert viewer["stains"] == ["A", "B"]
    plane = viewer["planes"][0]
    assert len(plane["fixed"]) == 2
    assert plane["fixed"][0] != plane["fixed"][1]
    assert len(plane["candidates"][0]["warped"]) == 1
    assert plane["candidates"][0]["warpedFor"] == [0, 0]
    assert len(viewer["stats"][0]["ncc"]) == 2


def test_report_hides_a_candidate_without_an_image_for_a_stain(tmp_path):
    fixed_a, fixed_b, moved_a, moved_b = _images(tmp_path)
    viewer = build_viewer_data(
        fixed_a,
        [{"name": "affine", "warped": moved_a}, {"name": "deform", "warped": moved_a}],
        [
            {"name": "A", "fixed": fixed_a, "warped": [moved_a, moved_a]},
            {"name": "B", "fixed": fixed_b, "warped": [None, moved_b]},
        ],
    )
    affine, deform = viewer["planes"][0]["candidates"]
    assert affine["warpedFor"] == [0, None]
    assert deform["warpedFor"] == [0, 1]
    assert viewer["stats"][0]["ncc"][1] is None
    assert viewer["stats"][1]["ncc"][1] is not None


def test_report_without_stains_has_one_unnamed_stain(tmp_path):
    fixed, _, moved, _ = _images(tmp_path)
    output = tmp_path / "report.html"
    render_report(
        str(output), "t", "e", {}, fixed, [{"name": "m", "warped": moved}]
    )
    assert "NCC to fixed" in output.read_text()
