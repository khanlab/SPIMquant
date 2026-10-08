"""The bottleneck conditioner: lantern_model.BottleneckConditioner and its wiring.

The conditioner feeds the network the two numbers per-tile normalisation
removes -- the log mean and log std of the RAW tile.  These tests cover the
arithmetic, the detection of whether a checkpoint wants it, and the fact that
the statistics reaching it are measured before normalisation (the bug that
would otherwise be silent: the model would read a standardised 0/1 for every
tile and the conditioning would carry no information).

Loads the real scripts, so it needs the runtime dependencies (torch,
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


def _load(name):
    """One of the workflow scripts, with its siblings importable while it loads."""
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        spec = spec_from_file_location(name, SCRIPTS_DIR / f"{name}.py")
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS_DIR))
    return module


@pytest.fixture(scope="module")
def lm():
    return _load("lantern_model")


@pytest.fixture(scope="module")
def lp():
    return _load("lantern_plaques")


# --------------------------------------------------------------- the arithmetic


def test_conditioner_standardises_the_log_statistics(lm):
    """The conditioning vector is (log1p(x) - mu) / sd for each statistic."""
    cond = lm.BottleneckConditioner(channels=8, cond_dim=4).eval()
    with torch.no_grad():
        cond.log_mean_mu.fill_(6.0)
        cond.log_mean_sd.fill_(2.0)
        cond.log_std_mu.fill_(5.0)
        cond.log_std_sd.fill_(4.0)

    v = cond.cond_vector(torch.tensor([[800.0]]), torch.tensor([[500.0]]))

    assert v.shape == (1, 2)
    np.testing.assert_allclose(
        v.detach().numpy(),
        [[(np.log1p(800.0) - 6.0) / 2.0, (np.log1p(500.0) - 5.0) / 4.0]],
        rtol=1e-6,
    )


def test_conditioner_is_identity_initialised(lm):
    """Identity on the features, zero on the conditioning: switching it on changes nothing."""
    cond = lm.BottleneckConditioner(channels=16, cond_dim=4).eval()
    feat = torch.randn(2, 16, 2, 2, 2)

    with torch.no_grad():
        out = cond(feat, torch.tensor([[800.0], [9.0]]), torch.tensor([[500.0], [3.0]]))

    torch.testing.assert_close(out, feat)


def test_conditioner_preserves_shape_and_reads_both_statistics(lm):
    cond = lm.BottleneckConditioner(channels=16, cond_dim=4).eval()
    with torch.no_grad():  # break the identity init so the conditioning can land
        cond.mix.weight.normal_(std=0.1)
    feat = torch.randn(1, 16, 2, 2, 2)
    one, two = torch.tensor([[800.0]]), torch.tensor([[500.0]])

    with torch.no_grad():
        base = cond(feat, one, two)
        other_mean = cond(feat, one * 10, two)
        other_std = cond(feat, one, two * 10)

    assert base.shape == feat.shape
    assert not torch.allclose(base, other_mean)
    assert not torch.allclose(base, other_std)


def test_conditioner_survives_a_nonpositive_mean(lm):
    """log1p is clamped at 0, so a background-subtracted store cannot produce a NaN."""
    cond = lm.BottleneckConditioner(channels=8, cond_dim=4).eval()

    v = cond.cond_vector(torch.tensor([[-5.0]]), torch.tensor([[0.0]]))

    assert torch.isfinite(v).all()


# ------------------------------------------------------------ the network wiring


def _small_net(lm, conditioner):
    """A ResEncLUNet small enough to run on CPU is not possible -- the encoder is fixed --
    so these tests use the real one but never call its forward."""
    return lm.ResEncLUNet(
        num_classes=2, num_input_channels=1, conditioner=conditioner, cond_dim=32
    )


def test_unconditioned_net_has_no_conditioner_and_needs_no_statistics(lm):
    net = _small_net(lm, conditioner=False)

    assert net.cond is None
    assert not any(k.startswith("cond.") for k in net.state_dict())


def test_conditioned_net_state_dict_matches_the_published_layout(lm):
    """The names and shapes the published checkpoints carry (apooladi/lantern-abeta)."""
    net = _small_net(lm, conditioner=True)
    sd = net.state_dict()

    assert sd["cond.proj.weight"].shape == (32, 2)
    assert sd["cond.proj.bias"].shape == (32,)
    assert sd["cond.mix.weight"].shape == (320, 352, 1, 1, 1)
    assert sd["cond.mix.bias"].shape == (320,)
    for name in ("log_mean_mu", "log_mean_sd", "log_std_mu", "log_std_sd"):
        assert sd[f"cond.{name}"].shape == (1,)


def test_decode_refuses_to_run_a_conditioned_fold_without_statistics(lm):
    """The model card is explicit that the statistics are required, not optional."""
    net = _small_net(lm, conditioner=True)

    with pytest.raises(ValueError, match="raw_mean and raw_std"):
        net.decode([torch.zeros(1, 320, 2, 2, 2)])


class _BottleneckOnlyDecoder(torch.nn.Module):
    """Stands in for the real decoder by handing back the (possibly conditioned)
    bottleneck, so a test can see what the decoder was given."""

    def forward(self, skips):
        return skips[-1]


def test_decode_does_not_write_through_the_shared_skips(lm):
    """Every fold decodes the SAME skip list from the one shared encoder pass, so
    conditioning fold 0 must not corrupt what fold 1 sees."""
    net = _small_net(lm, conditioner=True)
    with torch.no_grad():
        net.cond.mix.weight.normal_(std=0.1)
    bottleneck = torch.randn(1, 320, 2, 2, 2)
    skips = [bottleneck]
    net.decoder = _BottleneckOnlyDecoder()

    with torch.no_grad():
        out = net.decode(skips, torch.tensor([[800.0]]), torch.tensor([[500.0]]))

    assert not torch.allclose(out, bottleneck)  # the conditioner did something
    assert skips[0] is bottleneck  # ... to a copy


# ------------------------------------------------- reading it off the checkpoint


def _ckpt(cond_dim=None, n_cond=2, declare=None):
    """A checkpoint stub carrying only what checkpoint_cond_dim reads."""
    state = {"decoder.seg_layers.0.weight": torch.zeros(2, 320, 1, 1, 1)}
    if cond_dim is not None:
        state["cond.proj.weight"] = torch.zeros(cond_dim, n_cond)
    ckpt = {"state_dict": state, "num_classes": 2}
    declared = (cond_dim is not None) if declare is None else declare
    if declared:
        ckpt["bottleneck_cond"] = {"inputs": ["log1p(raw_mean)", "log1p(raw_std)"]}
    return ckpt


def test_cond_dim_is_none_for_an_unconditioned_checkpoint(lp):
    """The older ensembles (ki3-abeta, mapt_a/b/c-abeta) must keep loading."""
    assert lp.checkpoint_cond_dim(_ckpt(), "fold0") is None


def test_cond_dim_comes_from_the_weights(lp):
    assert lp.checkpoint_cond_dim(_ckpt(cond_dim=32), "fold0") == 32
    assert lp.checkpoint_cond_dim(_ckpt(cond_dim=8), "fold0") == 8


@pytest.mark.parametrize(
    "ckpt",
    [_ckpt(cond_dim=32, declare=False), _ckpt(declare=True)],
    ids=["weights but no metadata", "metadata but no weights"],
)
def test_cond_dim_rejects_a_checkpoint_that_contradicts_its_own_metadata(lp, ckpt):
    with pytest.raises(ValueError, match="bottleneck_cond"):
        lp.checkpoint_cond_dim(ckpt, "fold0")


def test_cond_dim_rejects_a_multi_stain_conditioner(lp):
    """Upstream's 4-slot conditioner takes 12 inputs; this pipeline feeds one stain."""
    with pytest.raises(ValueError, match="single stain"):
        lp.checkpoint_cond_dim(_ckpt(cond_dim=32, n_cond=12), "fold0")


# ------------------------------------------- the statistics reaching the network


def test_raw_tile_stats_are_the_raw_moments(lp):
    rng = np.random.default_rng(0)
    arr = (rng.random((3, 4, 4, 4), dtype=np.float32) * 400 + 700).astype(np.float32)

    mean, std = lp.raw_tile_stats(torch.from_numpy(arr))

    assert mean.shape == (3, 1) and std.shape == (3, 1)
    np.testing.assert_allclose(
        mean.numpy()[:, 0], arr.reshape(3, -1).mean(axis=1), rtol=1e-5
    )
    # ddof=0, numpy's convention and the one the checkpoints were standardised with
    np.testing.assert_allclose(
        std.numpy()[:, 0], arr.reshape(3, -1).std(axis=1), rtol=1e-5
    )


class _StatRecordingFold:
    """A fold that records the statistics it was handed and votes background."""

    n_classes = 2
    cond = object()

    def __init__(self):
        self.seen = []

    def encoder(self, x):
        return x

    def decode(self, skips, raw_mean=None, raw_std=None):
        self.seen.append((raw_mean.clone(), raw_std.clone()))
        return torch.cat(
            [torch.zeros_like(skips), torch.full_like(skips, -10.0)], dim=1
        )


def test_predict_volume_conditions_on_the_raw_not_the_normalised_tile(lp):
    """The whole point of the conditioner: it must see the tile's real brightness.

    A single tile covering the volume, deliberately far from zero mean and unit
    variance -- which is what the tile would look like after normalize_tiles.
    """
    fold = _StatRecordingFold()
    vol = (np.arange(64, dtype=np.float32).reshape(4, 4, 4) * 10.0) + 700.0

    lp.predict_volume(
        vol, [fold], torch.device("cpu"), tile=4, stride=4, batch_size=1,
        merge_mode="average",
    )

    assert len(fold.seen) == 1
    mean, std = fold.seen[0]
    assert mean.shape == (1, 1) and std.shape == (1, 1)
    np.testing.assert_allclose(float(mean), vol.mean(), rtol=1e-5)
    np.testing.assert_allclose(float(std), vol.std(), rtol=1e-5)


def test_predict_volume_measures_each_tile_separately(lp):
    """Two non-overlapping tiles of very different brightness get their own statistics."""
    fold = _StatRecordingFold()
    vol = np.empty((4, 4, 8), dtype=np.float32)
    vol[:, :, :4] = 100.0
    vol[:, :, 4:] = 5000.0

    lp.predict_volume(
        vol, [fold], torch.device("cpu"), tile=4, stride=4, batch_size=2,
        merge_mode="average",
    )

    means = sorted(
        float(v) for mean, _ in fold.seen for v in mean.flatten()
    )
    assert means == pytest.approx([100.0, 5000.0])
