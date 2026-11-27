"""The model-independent tiled-inference engine (tiled_inference.py).

Uses toy adapters, so no model weights are needed; it does need torch and
zarrnii, so run with `pixi run -e dev pytest`.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

torch = pytest.importorskip("torch")


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


SCRIPTS_DIR = _find_repo_root(Path(__file__).parent) / "spimquant/workflow/scripts"


@pytest.fixture(scope="module")
def ti():
    # Left on sys.path for the module: run_inference's spawned reader processes
    # import zarrnii_compat from it.
    sys.path.insert(0, str(SCRIPTS_DIR))
    import tiled_inference

    yield tiled_inference
    sys.path.remove(str(SCRIPTS_DIR))


def _identity(tiles):
    return tiles[None]


def test_config_from_overlap(ti):
    assert ti.TilingConfig.from_overlap(128, 0.5).stride == 64
    assert ti.TilingConfig.from_overlap(128, 0.0).stride == 128
    assert ti.TilingConfig.from_overlap(128, 0.25).halo == 32


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (dict(merge_mode="median"), "merge_mode"),
        (dict(stride=0), "stride"),
        (dict(tile=128, stride=256), "stride"),
        (dict(tile=128, chunk=64), "chunk"),
    ],
)
def test_config_rejects_bad_values(ti, kwargs, match):
    with pytest.raises(ValueError, match=match):
        ti.TilingConfig(**kwargs)


@pytest.mark.parametrize("mode", ["average", "gaussian", "max"])
def test_identity_model_reproduces_the_block(ti, mode):
    """Every merge mode returns the input when every tile returns its input."""
    vol = np.random.default_rng(0).random((10, 12, 9), dtype=np.float32)
    out = ti.predict_block(
        vol,
        _identity,
        n_out=1,
        device=torch.device("cpu"),
        tile=4,
        stride=2,
        batch_size=3,
        merge_mode=mode,
        tile_dtype=np.float32,
    )
    np.testing.assert_allclose(out[0], vol, rtol=1e-6)


def test_blocks_with_halo_match_one_whole_block(ti):
    """Splitting into halo-grown blocks changes nothing: no seams at block edges.

    The toy model is position-dependent (it scales each voxel by its distance
    from the tile centre), so a voxel's value depends on which tiles covered
    it; the halo must make each block see the same tiles as the whole volume.
    """
    tile, stride = 8, 4
    centre = np.abs(np.arange(tile) - (tile - 1) / 2)
    ramp = torch.from_numpy(
        (
            1 + centre[:, None, None] + centre[None, :, None] + centre[None, None, :]
        ).astype(np.float32)
    )

    def model(tiles):
        return (tiles * ramp)[None]

    shape = (1, 40, 36, 28)
    vol = np.random.default_rng(1).random(shape[1:], dtype=np.float32)
    kw = dict(
        n_out=1,
        device=torch.device("cpu"),
        tile=tile,
        stride=stride,
        batch_size=4,
        merge_mode="gaussian",
    )
    whole = ti.predict_block(vol, model, **kw)

    chunks = ti.plan_chunks(shape, 16, tile - stride)
    stitched = np.zeros_like(whole)
    plan = ti.block_plan(shape[1:], chunks[1:], tile - stride)
    assert len(plan) > 1
    for core, outer, inner in plan:
        part = ti.predict_block(vol[outer], model, **kw)
        stitched[(slice(None),) + core] = part[(slice(None),) + inner]
    np.testing.assert_allclose(stitched, whole, rtol=1e-5)


def test_plan_chunks_absorbs_a_tail_thinner_than_the_halo(ti):
    assert ti.plan_chunks((1, 1030, 512, 100), 512, 64) == (
        (1,),
        (512, 518),
        (512,),
        (100,),
    )


def test_device_pool_loads_once_per_device_and_lends_exclusively(ti):
    loaded = []
    pool = ti.DevicePool(["a", "b"], lambda d: loaded.append(d) or f"model-{d}")
    assert loaded == ["a", "b"]
    with pool.acquire() as (d1, m1), pool.acquire() as (d2, m2):
        assert {d1, d2} == {"a", "b"} and m1 == f"model-{d1}"


class _Doubler:
    """Minimal adapter: one channel, twice the input."""

    out_labels = ["double"]
    tile_dtype = np.float32
    n_out = 1

    def load(self, device):
        return 2.0

    def predict_tiles(self, model, tiles):
        return (tiles * model)[None]

    def finalize(self, merged):
        return merged

    def describe(self):
        return "doubler"


def test_run_inference_end_to_end(ti, tmp_path):
    """From an OME-Zarr on disk, through reader processes and blocks, to output."""
    import dask.array as da
    from zarrnii import ZarrNii

    img = np.random.default_rng(2).random((1, 48, 40, 36), dtype=np.float32)
    path = str(tmp_path / "in.ome.zarr")
    ZarrNii.from_darr(da.from_array(img), axes_order="ZYX").to_ome_zarr(path)
    spec = dict(path=path, level=0, target_scale=1000.0)
    from zarrnii_compat import open_inference_array

    data = open_inference_array(**spec)
    assert data.shape == img.shape

    config = ti.TilingConfig(tile=8, stride=4, batch_size=4, chunk=16, n_gpus=0)
    with ti.run_inference(_Doubler(), spec, config, threads=2, tmpdir=tmp_path) as out:
        result = np.asarray(out)
    np.testing.assert_allclose(result, 2 * np.asarray(data), rtol=1e-5)


def test_vesselfm_scaling_matches_monai(ti):
    """The batched torch rescale equals vesselFM's per-image MONAI transform."""
    pytest.importorskip("monai")
    from monai.transforms import ScaleIntensityRangePercentiles

    from importlib.util import module_from_spec, spec_from_file_location

    # Load the workflow script under another name: as "vesselfm" it would
    # shadow the vesselfm package it imports from.
    spec = spec_from_file_location("vesselfm_script", SCRIPTS_DIR / "vesselfm.py")
    script = module_from_spec(spec)
    spec.loader.exec_module(script)
    scale_percentiles = script.scale_percentiles

    rng = np.random.default_rng(3)
    tiles = rng.gamma(2.0, 300.0, size=(3, 16, 16, 16)).astype(np.float32)
    monai = ScaleIntensityRangePercentiles(
        lower=1, upper=99, b_min=0, b_max=1, clip=True
    )
    expected = np.stack([np.asarray(monai(t)) for t in tiles])
    got = scale_percentiles(torch.from_numpy(tiles)).numpy()
    np.testing.assert_allclose(got, expected, atol=1e-5)
