from contextlib import contextmanager
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import types

import numpy as np


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


def _load_lantern_module():
    repo_root = _find_repo_root(Path(__file__).parent)
    script_path = repo_root / "spimquant/workflow/scripts/lantern_plaques.py"

    dask = types.ModuleType("dask")
    dask_array = types.ModuleType("dask.array")
    dask_diagnostics = types.ModuleType("dask.diagnostics")

    class ProgressBar:  # pragma: no cover - only present to satisfy import
        pass

    dask_diagnostics.ProgressBar = ProgressBar
    dask.array = dask_array
    dask.diagnostics = dask_diagnostics
    sys.modules.setdefault("dask", dask)
    sys.modules.setdefault("dask.array", dask_array)
    sys.modules.setdefault("dask.diagnostics", dask_diagnostics)

    dask_setup = types.ModuleType("dask_setup")

    @contextmanager
    def get_dask_client(*args, **kwargs):
        yield None

    dask_setup.get_dask_client = get_dask_client
    sys.modules.setdefault("dask_setup", dask_setup)

    lantern_model = types.ModuleType("lantern_model")

    class ResEncLUNet:  # pragma: no cover - only present to satisfy import
        pass

    lantern_model.ResEncLUNet = ResEncLUNet
    sys.modules.setdefault("lantern_model", lantern_model)

    torch = types.ModuleType("torch")
    sys.modules.setdefault("torch", torch)

    zarrnii = types.ModuleType("zarrnii")

    class ZarrNii:  # pragma: no cover - only present to satisfy import
        pass

    zarrnii.ZarrNii = ZarrNii
    sys.modules.setdefault("zarrnii", zarrnii)

    spec = spec_from_file_location("lantern_plaques", script_path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lantern_plaques = _load_lantern_module()


def _merge_two_tiles(tile, volume_shape, first_votes, second_votes, merge_mode):
    votes = np.zeros(
        volume_shape,
        dtype=np.uint8 if merge_mode == "max" else np.float32,
    )
    weight_accum = (
        None if merge_mode == "max" else np.zeros(volume_shape, dtype=np.float32)
    )
    importance_map = (
        lantern_plaques.gaussian_importance_map(tile)
        if merge_mode == "gaussian"
        else None
    )
    batch_votes = np.stack([first_votes, second_votes])
    batch = [(0, 0, 0), (0, 0, 2)]
    lantern_plaques.merge_tile_votes(
        votes,
        batch_votes,
        batch,
        tile,
        merge_mode,
        n_folds=5,
        weight_accum=weight_accum,
        importance_map=importance_map,
    )
    if merge_mode == "max":
        return votes
    return votes / np.maximum(weight_accum, 1e-8)


def test_merge_mode_max_matches_existing_per_voxel_max_behavior():
    tile = 3
    volume_shape = (3, 3, 5)
    first_votes = np.full((tile, tile, tile), 1, dtype=np.uint8)
    second_votes = np.full((tile, tile, tile), 4, dtype=np.uint8)

    merged = _merge_two_tiles(tile, volume_shape, first_votes, second_votes, "max")

    expected = np.zeros(volume_shape, dtype=np.uint8)
    expected[:, :, :3] = 1
    expected[:, :, 2:5] = np.maximum(expected[:, :, 2:5], 4)
    np.testing.assert_array_equal(merged, expected)


def test_merge_mode_average_uniformly_averages_overlapping_tile_fractions():
    tile = 3
    volume_shape = (3, 3, 5)
    first_votes = np.full((tile, tile, tile), 5, dtype=np.uint8)
    second_votes = np.zeros((tile, tile, tile), dtype=np.uint8)

    merged = _merge_two_tiles(tile, volume_shape, first_votes, second_votes, "average")

    np.testing.assert_allclose(merged[:, :, 0:2], 1.0)
    np.testing.assert_allclose(merged[:, :, 2], 0.5)
    np.testing.assert_allclose(merged[:, :, 3:5], 0.0)


def test_merge_mode_gaussian_weights_tile_center_more_than_uniform_average():
    tile = 5
    volume_shape = (5, 5, 7)
    first_votes = np.full((tile, tile, tile), 5, dtype=np.uint8)
    second_votes = np.zeros((tile, tile, tile), dtype=np.uint8)

    average = _merge_two_tiles(tile, volume_shape, first_votes, second_votes, "average")
    gaussian = _merge_two_tiles(
        tile, volume_shape, first_votes, second_votes, "gaussian"
    )

    # Global x=2 lies at the center of the first tile and the edge of the second.
    assert gaussian[2, 2, 2] > average[2, 2, 2]
    # Global x=4 lies at the edge of the first tile and the center of the second.
    assert gaussian[2, 2, 4] < average[2, 2, 4]
