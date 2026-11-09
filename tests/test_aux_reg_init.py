"""Tests for the initial translation search of additional scans (aux_reg_init.py)."""

import importlib.util
import pathlib

import numpy as np
import pytest
from scipy import ndimage

_SCRIPT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "spimquant"
    / "workflow"
    / "scripts"
    / "aux_reg_init.py"
)
_spec = importlib.util.spec_from_file_location("aux_reg_init", _SCRIPT)
aux_reg_init = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(aux_reg_init)


def _blobs(shape, seed):
    rng = np.random.default_rng(seed)
    vol = ndimage.gaussian_filter(rng.random(shape).astype(np.float32), 2)
    return (vol - vol.min()) * 1000


def _affine(spacing, origin, flip_x=False):
    aff = np.diag([*spacing, 1.0])
    if flip_x:
        aff[0, 0] *= -1
    aff[:3, 3] = origin
    return aff


@pytest.mark.parametrize("flip_x", [False, True])
def test_match_finds_subvolume_with_unrelated_physical_origin(flip_x):
    """A finer, cropped scan with an arbitrary origin is located in the reference."""
    moving = _blobs((80, 70, 60), seed=0)
    moving_aff = _affine((0.1, 0.1, 0.1), (-3.0, 2.0, 5.0), flip_x)

    # the additional scan: the reference's voxels [20:44, 15:39, 10:30] at 2x
    # finer spacing, with an origin that has nothing to do with the reference's
    crop = moving[20:44, 15:39, 10:30]
    fixed = ndimage.zoom(crop, 2, order=1)
    true_corner = (moving_aff @ np.array([20, 15, 10, 1.0]))[:3]
    fixed_aff = _affine((0.05, 0.05, 0.05), (40.0, -12.0, 7.0), flip_x)
    # place the fine grid's first voxel where the crop's first voxel is, up to the
    # arbitrary origin offset (zoom keeps the corner voxel centres aligned)
    true_t = true_corner - fixed_aff[:3, 3]

    t = aux_reg_init.match_translation(
        [fixed], fixed_aff, [moving], moving_aff, min_res=0.1
    )
    np.testing.assert_allclose(t, true_t, atol=0.1 + 1e-6)


def test_center_translation_aligns_centers():
    fixed_aff = _affine((0.05,) * 3, (1.0, 2.0, 3.0))
    moving_aff = _affine((0.1,) * 3, (0.0, 0.0, 0.0))
    t = aux_reg_init.center_translation(fixed_aff, (21, 21, 21), moving_aff, (11, 11, 11))
    np.testing.assert_allclose(t, [-1.0, -2.0, -3.0])
