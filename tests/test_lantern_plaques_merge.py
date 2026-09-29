"""Tile-overlap merging (--merge_mode) in lantern_plaques.py.

Loads the real script, so it needs the runtime dependencies (torch,
dynamic-network-architectures); run with `pixi run -e dev pytest`.
"""

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("dynamic_network_architectures")


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


SCRIPTS_DIR = _find_repo_root(Path(__file__).parent) / "spimquant/workflow/scripts"


@pytest.fixture(scope="module")
def lp():
    """lantern_plaques.py, with its sibling helper modules importable while it loads."""
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        spec = spec_from_file_location(
            "lantern_plaques", SCRIPTS_DIR / "lantern_plaques.py"
        )
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS_DIR))
    return module


def _merge(lp, tile_votes, origins, shape, tile, mode):
    """merge_tile_votes over hand-made per-tile counts, normalised as predict_volume does.

    `tile_votes` is (n_out, n_tiles, tile, tile, tile); `origins` the tiles' corners.
    """
    n_out = tile_votes.shape[0]
    if mode == "max":
        votes = np.zeros((n_out,) + shape, dtype=np.uint8)
        weight = importance = None
    else:
        votes = np.zeros((n_out,) + shape, dtype=np.float32)
        weight = np.zeros(shape, dtype=np.float32)
        importance = lp.gaussian_importance_map(tile) if mode == "gaussian" else None
    lp.merge_tile_votes(votes, weight, tile_votes, origins, tile, importance)
    if weight is not None:
        np.divide(votes, weight, out=votes, where=weight > 0)
    return votes


def _two_tiles_along_x(first, second, tile=3, n_out=1):
    """Two tiles overlapping in the single x column 2, with constant vote counts."""
    tile_votes = np.zeros((n_out, 2, tile, tile, tile), dtype=np.uint8)
    tile_votes[:, 0] = np.asarray(first, dtype=np.uint8).reshape(n_out, 1, 1, 1)
    tile_votes[:, 1] = np.asarray(second, dtype=np.uint8).reshape(n_out, 1, 1, 1)
    origins = [(0, 0, 0), (0, 0, tile - 1)]
    shape = (tile, tile, 2 * tile - 1)
    return tile_votes, origins, shape


def test_max_keeps_the_highest_count(lp):
    tile_votes, origins, shape = _two_tiles_along_x(1, 4)

    merged = _merge(lp, tile_votes, origins, shape, tile=3, mode="max")

    assert merged.dtype == np.uint8
    np.testing.assert_array_equal(merged[0, :, :, :2], 1)
    np.testing.assert_array_equal(merged[0, :, :, 2:], 4)


def test_average_takes_the_mean_over_covering_tiles(lp):
    tile_votes, origins, shape = _two_tiles_along_x(5, 0)

    merged = _merge(lp, tile_votes, origins, shape, tile=3, mode="average")

    assert merged.dtype == np.float32
    np.testing.assert_allclose(merged[0, :, :, :2], 5.0)
    np.testing.assert_allclose(merged[0, :, :, 2], 2.5)
    np.testing.assert_allclose(merged[0, :, :, 3:], 0.0)


def test_gaussian_trusts_the_tile_whose_centre_is_nearer(lp):
    # tiles at x=0..4 and x=2..6, overlapping over x=2..4
    tile = 5
    tile_votes = np.zeros((1, 2, tile, tile, tile), dtype=np.uint8)
    tile_votes[:, 0] = 5
    origins = [(0, 0, 0), (0, 0, 2)]
    shape = (tile, tile, 7)

    average = _merge(lp, tile_votes, origins, shape, tile, mode="average")
    gaussian = _merge(lp, tile_votes, origins, shape, tile, mode="gaussian")

    # x=2 is the centre of the voting tile and the edge of the silent one ...
    assert gaussian[0, 2, 2, 2] > average[0, 2, 2, 2]
    # ... and x=4 the edge of the voting tile and the centre of the silent one.
    assert gaussian[0, 2, 2, 4] < average[0, 2, 2, 4]


def test_average_outvotes_a_single_edge_tile(lp):
    """The tiling artifact: one of the 8 tiles covering a voxel calls it 5/5."""
    tile, stride = 4, 2
    origins = [(z, y, x) for z in (0, 2) for y in (0, 2) for x in (0, 2)]
    tile_votes = np.zeros((1, len(origins), tile, tile, tile), dtype=np.uint8)
    tile_votes[0, 0] = 5
    shape = (tile + stride,) * 3
    centre = (0, 2, 2, 2)  # covered by all 8 tiles

    maxed = _merge(lp, tile_votes, origins, shape, tile, mode="max")
    averaged = _merge(lp, tile_votes, origins, shape, tile, mode="average")

    assert maxed[centre] == 5
    assert averaged[centre] == pytest.approx(5 / 8)
    # binarize_lantern_plaques cuts "3 of 5" at a mean of 2.5
    assert averaged[centre] < 2.5


def test_average_merges_every_channel_with_shared_weights(lp):
    # plaque channel: 4 and 2 votes; class-2 channel: 1 and 0 votes
    tile_votes, origins, shape = _two_tiles_along_x([4, 1], [2, 0], n_out=2)

    merged = _merge(lp, tile_votes, origins, shape, tile=3, mode="average")

    np.testing.assert_allclose(merged[0, :, :, 2], 3.0)
    np.testing.assert_allclose(merged[1, :, :, 2], 0.5)
    np.testing.assert_allclose(merged[1, :, :, :2], 1.0)
    np.testing.assert_allclose(merged[1, :, :, 3:], 0.0)


class _FakeFold:
    """Stands in for one ResEncLUNet fold that calls every voxel of every tile
    the same way: plaque if `votes_plaque`, background otherwise."""

    n_classes = 2

    def __init__(self, votes_plaque):
        self.logit = 10.0 if votes_plaque else -10.0

    def encoder(self, x):
        return x

    def decoder(self, skips):
        background = torch.zeros_like(skips)
        return torch.cat([background, torch.full_like(skips, self.logit)], dim=1)


@pytest.mark.parametrize("mode", ["max", "average", "gaussian"])
def test_predict_volume_ignores_tiles_skipped_for_no_brain(lp, mode):
    nets = [_FakeFold(True)] * 3 + [_FakeFold(False)] * 2
    vol = np.random.default_rng(0).random((4, 4, 8), dtype=np.float32)
    brain = np.zeros(vol.shape, dtype=bool)
    brain[:, :, :2] = True  # only the tile at x=0 (of x=0, 2, 4) touches brain

    votes = lp.predict_volume(
        vol,
        nets,
        torch.device("cpu"),
        tile=4,
        stride=2,
        batch_size=2,
        brain=brain,
        merge_mode=mode,
    )

    assert votes.shape == (1,) + vol.shape
    assert votes.dtype == (np.uint8 if mode == "max" else np.float32)
    # the skipped tiles neither vote nor dilute the average where they overlap
    np.testing.assert_allclose(votes[0, :, :, :4], 3.0, rtol=1e-6)
    np.testing.assert_array_equal(votes[0, :, :, 4:], 0)


def test_predict_volume_rejects_an_unknown_merge_mode(lp):
    with pytest.raises(ValueError, match="merge_mode"):
        lp.predict_volume(
            np.zeros((4, 4, 4), dtype=np.float32),
            [_FakeFold(True)],
            torch.device("cpu"),
            tile=4,
            stride=2,
            batch_size=1,
            merge_mode="median",
        )
