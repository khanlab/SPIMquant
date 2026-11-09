"""End-to-end test of the fibre_moments and fibre_metrics workflow scripts."""

import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

zarrnii = pytest.importorskip("zarrnii")
nib = pytest.importorskip("nibabel")


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


ROOT = _find_repo_root(Path(__file__).parent)
SCRIPTS_DIR = ROOT / "spimquant/workflow/scripts"


def _load(name):
    spec = spec_from_file_location(name, SCRIPTS_DIR / f"{name}.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def scripts():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        yield _load("fibre_moments"), _load("fibre_metrics")
    finally:
        sys.path.remove(str(SCRIPTS_DIR))


def _phantom(shape, spacing_um, seed=0):
    """Fibres along array axis 2 in one half, along axis 1 in the other,
    inside an ellipsoidal 'brain' (zero outside, as the corrected image)."""
    rng = np.random.default_rng(seed)
    grids = np.meshgrid(
        *[np.arange(n) * s for n, s in zip(shape, spacing_um)], indexing="ij"
    )
    coords = np.stack(grids, axis=-1)
    ext = np.asarray(shape) * np.asarray(spacing_um)
    img = np.zeros(shape, np.float32)
    for k in range(60):
        axis = 2 if k % 2 == 0 else 1
        d = np.zeros(3)
        d[axis] = 1.0
        p = rng.uniform(0.1, 0.9, 3) * ext
        r = coords - p
        perp = r - (r @ d)[..., None] * d
        line = 100 * np.exp(-np.sum(perp**2, -1) / (2 * 6.0**2))
        half = (
            coords[..., 0] < ext[0] / 2 if axis == 2 else coords[..., 0] >= ext[0] / 2
        )
        img = np.maximum(img, line * half)
    img += 20 + rng.normal(0, 2, shape)
    centre = ext / 2
    inside = np.sum(((coords - centre) / (0.48 * ext)) ** 2, -1) <= 1
    return np.where(inside, img, 0).astype(np.float32)


def test_moments_and_metrics(scripts, tmp_path):
    import dask.array as da

    moments_mod, metrics_mod = scripts
    spacing_um = (4.0, 4.0, 4.0)
    shape = (64, 64, 64)
    img = _phantom(shape, spacing_um)
    corrected = tmp_path / "corrected.ome.zarr"
    zarrnii.ZarrNii.from_darr(
        da.from_array(img[None], chunks=(1, 32, 32, 32)),
        spacing=tuple(s / 1000 for s in spacing_um),
        origin=(1.0, 2.0, 3.0),
        axes_units={"x": "millimeter", "y": "millimeter", "z": "millimeter"},
    ).to_ome_zarr(str(corrected), max_layer=3)

    out = SimpleNamespace(
        moments=str(tmp_path / "moments.ome.zarr"),
        params_json=str(tmp_path / "params.json"),
        qc_png=str(tmp_path / "qc.png"),
    )
    moments_mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(corrected=str(corrected)),
        output=out,
        params=SimpleNamespace(
            scales_um=[6, 10],
            block_um=32,
            lmax=4,
            alpha=0.001,
            beta=1.0,
            gamma_quantile=0.999,
            gamma_factor=0.5,
            threshold_method="yen",
            stat_blocks=4,
            stat_samples_per_block=5000,
            stat_level_offset=2,
            edge_min_weight=0.5,
            chunk_size=32,
            dask_threads_per_worker=1,
            dask_reserve_mb=0,
        ),
        threads=2,
        resources={},
        config={"dask_scheduler": "threads", "zarrnii_out_kwargs": {}},
    )
    moments_mod.main()

    params = json.loads(Path(out.params_json).read_text())
    assert params["block_factor"] == [8, 8, 8]
    assert len(params["gammas"]) == 2
    assert Path(out.qc_png).exists()

    names = list(metrics_mod.SCALAR_METRICS) + ["dir", "dec", "odf"]
    mout = SimpleNamespace(**{n: str(tmp_path / f"{n}.nii.gz") for n in names})
    metrics_mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(moments=out.moments), output=mout, threads=2
    )
    metrics_mod.main()

    ref = nib.load(mout.fa)
    assert ref.shape == (8, 8, 8)
    # block grid: 8x the voxel size (mm), centred on the blocks
    np.testing.assert_allclose(np.abs(np.diag(ref.affine)[:3]), 0.032, rtol=1e-5)
    dens = np.asarray(nib.load(mout.fibredensity).dataobj)
    fa = np.asarray(ref.dataobj)
    assert dens.max() > 0.01
    # NIfTI is x, y, z: array axis 0 (zarr z) is NIfTI axis 2. Fibres along
    # array axis 2 (world x) in the lower half of axis 0, along axis 1
    # (world y) in the upper half.
    direction = np.asarray(nib.load(mout.dir).dataobj)
    assert direction.shape == (8, 8, 8, 3)
    has = (dens > 0.02) & (fa > 0.5)
    lower = has.copy()
    lower[:, :, 4:] = False
    upper = has.copy()
    upper[:, :, :4] = False
    assert lower.sum() > 5 and upper.sum() > 5
    assert np.median(np.abs(direction[..., 0][lower])) > 0.9
    assert np.median(np.abs(direction[..., 1][upper])) > 0.9
    odf = nib.load(mout.odf)
    assert odf.shape == (8, 8, 8, 15)
