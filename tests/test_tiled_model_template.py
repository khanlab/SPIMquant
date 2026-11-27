"""The new-model template (docs/howto/templates/tiled_model.py) still runs.

Fills the template's one TODO that has no default -- build_network -- with a
tiny 3D conv, saves it as a checkpoint, and runs the adapter through the
engine, so the template cannot silently drift from tiled_inference's API.
Needs torch and zarrnii: run with `pixi run -e dev pytest`.
"""

import sys
from importlib.util import module_from_spec, spec_from_file_location
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


ROOT = _find_repo_root(Path(__file__).parent)
SCRIPTS_DIR = ROOT / "spimquant/workflow/scripts"
TEMPLATE = ROOT / "docs/howto/templates/tiled_model.py"


@pytest.fixture(scope="module")
def tmpl():
    # Appended, not prepended, so that workflow scripts cannot shadow the
    # packages they import; left in place for the spawned reader processes.
    sys.path.append(str(SCRIPTS_DIR))
    spec = spec_from_file_location("tiled_model_template", TEMPLATE)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    yield module
    sys.path.remove(str(SCRIPTS_DIR))


def _tiny_net():
    torch.manual_seed(0)
    return torch.nn.Conv3d(1, 1, kernel_size=3, padding=1)


@pytest.fixture
def adapter(tmpl, tmp_path):
    class Tiny(tmpl.MyModel):
        def build_network(self):
            return _tiny_net()

    path = tmp_path / "tiny.pt"
    torch.save(_tiny_net().state_dict(), path)
    return Tiny(str(path), threshold=0.5)


def test_template_build_network_is_left_to_fill_in(tmpl):
    with pytest.raises(NotImplementedError):
        tmpl.MyModel("unused.pt").build_network()


def test_template_adapter_shapes_and_output(adapter):
    import dask.array as da

    net = adapter.load(torch.device("cpu"))
    tiles = torch.rand(3, 8, 8, 8)
    out = adapter.predict_tiles(net, tiles)
    assert tuple(out.shape) == (adapter.n_out, 3, 8, 8, 8)
    assert float(out.min()) >= 0.0 and float(out.max()) <= 1.0  # probabilities

    stored = adapter.finalize(da.from_array(np.array([[[[0.2, 0.7]]]], np.float32)))
    assert stored.dtype == np.uint8
    assert stored.compute().ravel().tolist() == [0, 100]


def test_template_runs_end_to_end_through_the_engine(tmpl, adapter, tmp_path):
    """From an OME-Zarr on disk, through reader processes, to the merged mask."""
    import dask.array as da
    from zarrnii import ZarrNii

    import tiled_inference as ti

    img = np.random.default_rng(0).random((1, 40, 36, 32), dtype=np.float32)
    path = str(tmp_path / "in.ome.zarr")
    ZarrNii.from_darr(da.from_array(img), axes_order="ZYX").to_ome_zarr(path)
    spec = dict(path=path, level=0, target_scale=None)

    config = ti.TilingConfig(tile=8, stride=4, batch_size=4, chunk=16, n_gpus=0)
    with ti.run_inference(adapter, spec, config, threads=2, tmpdir=tmp_path) as out:
        mask = np.asarray(out)
    assert mask.shape == (1,) + img.shape[1:]
    assert set(np.unique(mask)) <= {0, 100}
