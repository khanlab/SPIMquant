"""DoG tuning report helpers (dog_tune_report.py)."""

import sys
from pathlib import Path

import numpy as np
import pytest


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


SCRIPTS_DIR = _find_repo_root(Path(__file__).parent) / "spimquant/workflow/scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import dog_tune_report as dtr  # noqa: E402
import dogthresh as dt  # noqa: E402

ZOOMS = (5.0, 3.5, 3.5)


def test_weighted_dog_matches_dog_response_and_is_linear():
    rng = np.random.default_rng(0)
    img = rng.uniform(50, 150, (24, 40, 40)).astype(np.float32)
    img[:, :, :8] = 0
    s1, s2 = dt.sigma_vox(3, ZOOMS), dt.sigma_vox(12, ZOOMS)
    w = (img > 0).astype(np.float32)
    ref = dt.dog_response(img, s1, s2)[0]
    got = dtr.weighted_dog(img, w, s1, s2)[0]
    np.testing.assert_allclose(got, ref, atol=1e-4)
    other = rng.uniform(0, 10, img.shape).astype(np.float32) * w
    both = dtr.weighted_dog(img + 2.5 * other, w, s1, s2)[0]
    sep = got + 2.5 * dtr.weighted_dog(other, w, s1, s2)[0]
    np.testing.assert_allclose(both, sep, atol=1e-3)


def test_noise_model_recovers_white_noise_sd():
    rng = np.random.default_rng(1)
    img = (100 + 7 * rng.standard_normal((40, 64, 64))).astype(np.float32)
    tissue = np.ones(img.shape, bool)
    _, centres, sds = dtr.noise_model(img, tissue, ZOOMS)
    assert centres is not None
    np.testing.assert_allclose(sds, 7, rtol=0.1)


def test_noise_null_is_structureless():
    rng = np.random.default_rng(2)
    img = np.full((30, 60, 60), 100, np.float32)
    img[10:20, 20:40, 20:40] = 400  # bright structure
    img += rng.standard_normal(img.shape).astype(np.float32) * 5
    tissue = img > 0
    null, sd = dtr.noise_null_image(img, tissue, ZOOMS, rng)
    # constant mean everywhere: the bright square is gone
    assert abs(null[15, 30, 30] - null[2, 5, 5]) < 6 * 5
    assert sd == pytest.approx(5, rel=0.2)


def test_contrast_and_k_helpers_are_consistent():
    rng = np.random.default_rng(3)
    base = rng.normal(0, 1, 2000)
    gain = rng.uniform(4, 8, 2000)
    k = 6.0
    c90 = dtr.contrast_at_recall(base, gain, k, 0.9)
    assert dtr.recall_at(base, gain, k, c90 * 1.001) >= 0.9
    assert dtr.recall_at(base, gain, k, c90 * 0.98) < 0.9
    k90 = dtr.k_for_target(base, gain, 1.0, 0.9)
    assert dtr.recall_at(base, gain, k90 - 1e-6, 1.0) >= 0.9
    # never-detectable puncta push the limit to infinity
    assert dtr.contrast_at_recall(base, np.zeros_like(gain), k, 0.5) == np.inf


def test_spikes_are_detected_at_known_contrast():
    """A punctum at 100% contrast on flat tissue gives a DoG peak in line with
    the analytic Gaussian convolution."""
    shape = (40, 80, 80)
    w = np.ones(shape, np.float32)
    pts = np.array([[20, 40, 40]])
    sd_um = 4.0
    spk = dtr.spike_image(shape, pts, [100.0], sd_um, ZOOMS)
    assert spk[20, 40, 40] == pytest.approx(100)
    s = 3.0
    d = dtr.weighted_dog(spk, w, dt.sigma_vox(s, ZOOMS), dt.sigma_vox(4 * s, ZOOMS))[0]
    expect = 100 * (
        (sd_um**2 / (sd_um**2 + s**2)) ** 1.5
        - (sd_um**2 / (sd_um**2 + (4 * s) ** 2)) ** 1.5
    )
    assert d[20, 40, 40] == pytest.approx(expect, rel=0.1)


def test_place_spikes_respects_separation():
    rng = np.random.default_rng(4)
    allowed = np.ones((30, 60, 60), bool)
    pts = dtr.place_spikes(allowed, ZOOMS, 20, 40.0, rng)
    assert len(pts) > 3
    phys = pts * np.asarray(ZOOMS)
    d = np.linalg.norm(phys[:, None] - phys[None], axis=-1)
    assert d[np.triu_indices(len(pts), 1)].min() >= 40 - np.max(ZOOMS)


def test_scale_census_finds_blob_size():
    shape = (48, 96, 96)
    pts = np.array([[24, 30, 30], [24, 66, 66]])
    img = 50 + dtr.spike_image(shape, pts, [200.0, 200.0], 8.0, ZOOMS)
    valid = np.ones(shape, bool)
    blob_sd = dtr.scale_census(img, valid, ZOOMS, [2, 3, 4, 5, 6, 7, 8.5, 10, 12], 2)
    assert len(blob_sd) == 2
    np.testing.assert_allclose(blob_sd, 8.0, rtol=0.2)


def test_method_string():
    assert dtr.method_string(3.0, 4.0, 6.5) == "dog+s3r4k6p5"
    assert dt.parse_dog_method(dtr.method_string(2.5, 4, 6)) == (2.5, 4.0, 6.0)
