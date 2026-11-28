"""The segmentation overview QC reorients volumes to RAS before slicing.

Imaris and OME-Zarr inputs are stored ZYX and in the scan's own orientation;
the figure's axial/coronal/sagittal rows are only right if the volume is
reoriented with its affine first.
"""

import importlib.util
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "spimquant"
    / "workflow"
    / "scripts"
    / "qc_segmentation_overview.py"
)


@pytest.fixture(scope="module")
def overview():
    spec = importlib.util.spec_from_file_location("qc_segmentation_overview", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["qc_segmentation_overview"] = module
    spec.loader.exec_module(module)
    return module


def _marked_volume():
    """10 x 20 x 30 volume (x, y, z) with a bright voxel at the R, A, S corner."""
    data = np.zeros((10, 20, 30), dtype=np.float32)
    data[-1, -1, -1] = 1.0
    return data


def test_ras_volume_reorients_lpi(overview, tmp_path):
    from zarrnii import ZarrNii

    # Store the volume LPI: flip x and y, so the marked corner sits at index 0
    data = _marked_volume()[::-1, ::-1, :].copy()
    affine = np.diag([-0.1, -0.2, 0.3, 1.0])  # x->L, y->P, z->S, mm
    path = tmp_path / "lpi.nii.gz"
    nib.save(nib.Nifti1Image(data, affine), path)

    vol, zooms = overview._ras_volume(ZarrNii.from_nifti(str(path)))

    assert vol.shape == (10, 20, 30)
    assert np.allclose(zooms, [0.1, 0.2, 0.3])
    # back in RAS, the marked voxel is at the far R, A, S corner again
    assert np.unravel_index(np.argmax(vol), vol.shape) == (9, 19, 29)


@pytest.mark.parametrize(
    "step_axis, expected",
    [(2, 0.2 / 0.1), (1, 0.3 / 0.1), (0, 0.3 / 0.2)],
)
def test_slice_aspect(overview, step_axis, expected):
    assert overview._slice_aspect((0.1, 0.2, 0.3), step_axis) == pytest.approx(expected)
