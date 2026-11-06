"""DoG blob segmentation with pool exclusion (dogthresh.py)."""

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage as ndi


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


SCRIPTS_DIR = _find_repo_root(Path(__file__).parent) / "spimquant/workflow/scripts"


@pytest.fixture(scope="module")
def dt():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        spec = spec_from_file_location("dogthresh", SCRIPTS_DIR / "dogthresh.py")
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS_DIR))
    return module


# Voxel size of the marmoset data the defaults were tuned on (z, y, x), in um
ZOOMS = (5.0, 3.5, 3.5)


@pytest.mark.parametrize(
    "method, expected",
    [
        ("dog+s3r4k6", (3.0, 4.0, 6.0)),
        ("dog+s2p5r4k6p5", (2.5, 4.0, 6.5)),
        ("dog+s10r3p5k12", (10.0, 3.5, 12.0)),
    ],
)
def test_parse_dog_method(dt, method, expected):
    assert dt.parse_dog_method(method) == expected


@pytest.mark.parametrize(
    "method", ["dog+s3k6", "dog+s3r4", "gmm+n3k1", "dog+s3r1k6", "dog+s0r4k6"]
)
def test_parse_dog_method_rejects(dt, method):
    with pytest.raises(ValueError):
        dt.parse_dog_method(method)


def test_robust_threshold(dt):
    rng = np.random.default_rng(0)
    values = np.concatenate([rng.normal(0, 2, 100_000), np.full(1000, 500.0)])
    thr, med, sigma = dt.robust_threshold(values, 6)
    # MAD ignores the 1% of outliers
    assert abs(med) < 0.1
    assert sigma == pytest.approx(2.0, rel=0.05)
    assert thr == pytest.approx(med + 6 * sigma)


def _synthetic(rng):
    """Tissue at ~90 with noise, a few puncta and a large bright surface slab."""
    shape = (40, 160, 160)
    img = rng.normal(90, 5, shape).astype(np.float32)
    zz, yy, xx = np.indices(shape)
    centres = [(20, 40, 40), (20, 40, 100), (12, 110, 60), (28, 100, 110)]
    for cz, cy, cx in centres:
        r2 = ((zz - cz) * ZOOMS[0]) ** 2 + ((yy - cy) * ZOOMS[1]) ** 2
        r2 = r2 + ((xx - cx) * ZOOMS[2]) ** 2
        img += 250 * np.exp(-r2 / (2 * 5.0**2))  # ~12 um FWHM, dim puncta
    # one very bright punctum, as bright as the pool, but small
    r2 = ((zz - 20) * 5.0) ** 2 + ((yy - 80) * 3.5) ** 2 + ((xx - 40) * 3.5) ** 2
    img += 5000 * np.exp(-r2 / (2 * 5.0**2))
    # surface pool: a 25-voxel-thick slab at the far x edge, at 5000
    img[:, :, 130:155] = 5000 + rng.normal(0, 50, (shape[0], shape[1], 25))
    # outside the brain mask the corrected image is zero
    img[:, :, 155:] = 0
    return img, centres


def _run(dt, img, sigma_um=3.0, ratio=4.0, k=6.0):
    sigma1 = dt.sigma_vox(sigma_um, ZOOMS)
    sigma2 = dt.sigma_vox(sigma_um * ratio, ZOOMS)
    dog, _, tissue = dt.dog_response(img, sigma1, sigma2)
    thr, _, _ = dt.robust_threshold(dog[tissue], k)

    # pool mask on a 4x coarser grid, as the workflow does at level +2
    lo = img.reshape(10, 4, 40, 4, 40, 4).mean(axis=(1, 3, 5))
    lo_zooms = tuple(4 * z for z in ZOOMS)
    tissue_median = float(np.median(lo[lo > 0]))
    pool_thresh = 4 * tissue_median
    pool_lo = dt.pool_mask_lowres(lo, pool_thresh, 1e6, lo_zooms, 30)
    pool_up = dt.upsample_nearest_block(pool_lo, img.shape, (0, 0, 0), img.shape)

    mask = dt.segment_block(
        img,
        pool_up,
        sigma1,
        sigma2,
        thr,
        pool_thresh,
        dt.ellipsoid(7, ZOOMS),
    )
    return mask, pool_lo


def test_puncta_kept_pool_dropped(dt):
    img, centres = _synthetic(np.random.default_rng(1))
    mask, pool_lo = _run(dt, img)
    assert mask.dtype == np.uint8
    assert set(np.unique(mask)) <= {0, 100}

    for c in centres:
        assert mask[c] == 100, f"dim punctum at {c} missed"
    assert mask[20, 80, 40] == 100, "bright punctum dropped"

    assert pool_lo.any()
    # nothing in or along the pool
    assert not mask[:, :, 125:].any()

    # only puncta-sized objects remain
    labels, n = ndi.label(mask > 0)
    assert n == len(centres) + 1
    sizes = np.bincount(labels.ravel())[1:] * np.prod(ZOOMS)
    assert sizes.max() < 2e4


def test_no_rim_at_brain_mask_edge(dt):
    """The step from tissue to the zeroed outside must not respond."""
    rng = np.random.default_rng(2)
    img = rng.normal(90, 5, (30, 80, 80)).astype(np.float32)
    img[:, :, 50:] = 0
    sigma1 = dt.sigma_vox(3.0, ZOOMS)
    sigma2 = dt.sigma_vox(12.0, ZOOMS)
    dog, _, tissue = dt.dog_response(img, sigma1, sigma2)
    thr, _, _ = dt.robust_threshold(dog[tissue], 6)
    assert not (dog[:, :, 40:50] > thr).any()
    assert not dog[~tissue].any()


def test_upsample_nearest_block_matches_full(dt):
    lo = np.arange(5 * 6 * 7).reshape(5, 6, 7)
    full_shape = (21, 25, 29)  # not exact multiples, as in real pyramids
    full = dt.upsample_nearest_block(lo, full_shape, (0, 0, 0), full_shape)
    block = dt.upsample_nearest_block(lo, (8, 9, 10), (6, 7, 11), full_shape)
    np.testing.assert_array_equal(block, full[6:14, 7:16, 11:21])
    # corners map to the corresponding lowres corners
    assert full[0, 0, 0] == lo[0, 0, 0]
    assert full[-1, -1, -1] == lo[-1, -1, -1]


def test_select_stat_blocks(dt):
    fractions = {(0, 0, i): f for i, f in enumerate([0.0, 0.1, 0.6, 0.9, 1.0])}
    chosen = dt.select_stat_blocks(fractions, 2)
    assert set(chosen) <= {(0, 0, 2), (0, 0, 3), (0, 0, 4)}
    assert len(chosen) == 2
    # falls back to tissue-poor blocks, never to empty ones
    chosen = dt.select_stat_blocks(fractions, 10)
    assert set(chosen) == {(0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4)}


def test_tissue_fraction_per_block(dt):
    lo = np.zeros((4, 4, 4))
    lo[:, :, :2] = 1
    chunks = ((8, 8), (8, 8), (8, 8))
    fractions = dt.tissue_fraction_per_block(lo, chunks, (16, 16, 16))
    assert fractions[(0, 0, 0)] == 1.0
    assert fractions[(1, 1, 1)] == 0.0
