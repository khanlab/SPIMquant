"""Brain-masked intensities without bias-field correction (mask_correction.py)."""

import sys
from pathlib import Path

import dask.array as da
import numpy as np
import pytest


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


SCRIPTS_DIR = _find_repo_root(Path(__file__).parent) / "spimquant/workflow/scripts"


@pytest.fixture(scope="module")
def mc():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        import mask_correction
    finally:
        sys.path.remove(str(SCRIPTS_DIR))
    return mask_correction


@pytest.fixture(scope="module")
def dt():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        import dogthresh
    finally:
        sys.path.remove(str(SCRIPTS_DIR))
    return dogthresh


def _run(mc, img, lowres_mask, method, chunks):
    arr = da.from_array(img, chunks=chunks)
    out_dtype = np.float32 if method == "log" else img.dtype
    return da.map_blocks(
        mc.mask_block,
        arr,
        lowres_mask=lowres_mask,
        full_shape=img.shape,
        method=method,
        dtype=out_dtype,
        meta=np.array((), dtype=out_dtype),
    ).compute(scheduler="synchronous")


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    img = rng.integers(100, 4000, (32, 64, 64)).astype(np.uint16)
    lowres_mask = np.zeros((4, 8, 8), dtype=bool)  # 8x downsampled
    lowres_mask[1:3, 2:6, 2:6] = True
    return img, lowres_mask


@pytest.mark.parametrize("method", ["none", "log"])
def test_mask_block(mc, data, method):
    img, lowres_mask = data
    inside = mc.upsample_mask_block(lowres_mask, img.shape, (0, 0, 0), img.shape)
    # the box's interior is inside; far from it is outside
    assert inside[12:20, 20:44, 20:44].all()
    assert not inside[:4].any() and not inside[:, :, :12].any()
    out = _run(mc, img, lowres_mask, method, chunks=(16, 32, 32))
    assert not out[~inside].any()
    if method == "none":
        assert out.dtype == np.uint16
        np.testing.assert_array_equal(out[inside], img[inside])
    else:
        assert out.dtype == np.float32
        np.testing.assert_allclose(
            out[inside], np.log1p(img[inside].astype(float)), rtol=1e-6
        )
        # strictly positive inside, so downstream "img > 0 is tissue" still holds
        assert (out[inside] > 0).all()


def test_mask_block_independent_of_chunks(mc, data):
    img, lowres_mask = data
    a = _run(mc, img, lowres_mask, "none", chunks=(16, 32, 32))
    b = _run(mc, img, lowres_mask, "none", chunks=(7, 13, 50))
    np.testing.assert_array_equal(a, b)


def test_upsample_mask_block_is_smooth(mc):
    """A diagonal low-res edge becomes a fine staircase, not 8-voxel steps."""
    low = np.zeros((1, 8, 8), dtype=bool)
    low[0][np.tril_indices(8)] = True
    full = mc.upsample_mask_block(low, (1, 64, 64), (0, 0, 0), (1, 64, 64))[0]
    # boundary column per row moves by at most a few voxels between rows
    edge = full.sum(axis=1)
    assert np.abs(np.diff(edge[4:56].astype(int))).max() <= 2
    # and blocks agree with the whole
    part = mc.upsample_mask_block(low, (1, 20, 30), (0, 10, 17), (1, 64, 64))[0]
    np.testing.assert_array_equal(part, full[10:30, 17:47])


def test_check_mask_shape(mc):
    mc.check_mask_shape((30, 136, 115), (1909, 8731, 7386), 64)
    with pytest.raises(ValueError):
        mc.check_mask_shape((30, 136, 115), (1909, 8731, 7386), 32)


def test_pool_threshold(dt):
    assert dt.pool_threshold(90.0, 4) == 360.0
    # in log space, "4x the median" is the same ratio as an offset
    med = np.log1p(90.0)
    thr = dt.pool_threshold(med, 4, log_intensity=True)
    assert np.expm1(thr) == pytest.approx(4 * 90, rel=0.02)
