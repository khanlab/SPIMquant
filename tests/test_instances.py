"""Instance splitting of segmentation masks (instances.py)."""

import sys
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
def ins():
    # instances imports dogthresh from the scripts directory, as Snakemake
    # scripts do
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        import instances
    finally:
        sys.path.remove(str(SCRIPTS_DIR))
    return instances


# Voxel size of the marmoset data the defaults were tuned on (z, y, x), in um
ZOOMS = (5.0, 3.545, 3.545)


def _puncta(shape, centres_um, sigma_um=4.0, amp=200.0, noise=5.0, seed=0):
    rng = np.random.default_rng(seed)
    img = rng.normal(90, noise, shape).astype(np.float32)
    axes = [np.arange(n) * z for n, z in zip(shape, ZOOMS)]
    for c in centres_um:
        g = [np.exp(-((a - p) ** 2) / (2 * sigma_um**2)) for a, p in zip(axes, c)]
        img += amp * g[0][:, None, None] * g[1][None, :, None] * g[2][None, None, :]
    return img


def _mask(img, thresh=150):
    return (ndi.gaussian_filter(img, [1.0, 1.0, 1.0]) > thresh).astype(np.uint8) * 100


def test_neighbourhood_includes_face_neighbours(ins):
    fp = ins.neighbourhood(4.0, ZOOMS)  # 0.8 voxel in z
    assert fp.shape[0] == 3
    centre = tuple(s // 2 for s in fp.shape)
    for axis in range(3):
        for step in (-1, 1):
            idx = list(centre)
            idx[axis] += step
            assert fp[tuple(idx)]


def test_touching_puncta_split(ins):
    """Two puncta 5 sigma apart form one connected component but two instances."""
    shape = (24, 64, 64)
    ctr = np.array(shape) * ZOOMS / 2
    d = np.array([0, 0, 8.0])  # 16 um (4 sigma) apart along x
    img = _puncta(shape, [ctr - d, ctr + d])
    mask = _mask(img, thresh=120)
    assert ndi.label(mask)[1] == 1
    seeds, labels, _ = ins.segment_instances(img, mask, ZOOMS, 3.0, 4.0, h=7.0)
    assert len(seeds) == 2
    assert (labels[mask > 0] > 0).all()


def test_single_punctum_not_split(ins):
    shape = (24, 64, 64)
    img = _puncta(shape, [np.array(shape) * ZOOMS / 2], sigma_um=6.0)
    seeds, _, _ = ins.segment_instances(img, _mask(img), ZOOMS, 3.0, 4.0, h=7.0)
    assert len(seeds) == 1


def test_every_component_gets_a_seed(ins):
    """A faint, flat object below h still becomes one instance."""
    mask = np.zeros((10, 30, 30), dtype=np.uint8)
    mask[3:6, 5:12, 5:12] = 100
    mask[3:6, 18:25, 18:25] = 100
    img = np.full(mask.shape, 90, dtype=np.float32)
    img[mask > 0] = 95
    seeds, labels, _ = ins.segment_instances(img, mask, ZOOMS, 3.0, 4.0, h=50.0)
    assert len(seeds) == 2
    assert (labels[mask > 0] > 0).all()


def test_instance_features(ins):
    labels = np.zeros((10, 40, 40), dtype=np.int32)
    labels[4:6, 10:12, 5:35] = 1  # rod along x
    labels[2:5, 20:23, 20:23] = 2
    intensity = np.where(labels == 2, 10.0, 1.0)
    f = ins.instance_features(labels, [1, 2], intensity, ZOOMS)
    assert f["nvoxels"].tolist() == [120, 27]
    np.testing.assert_allclose(f["centroid"][1], [3, 21, 21])
    assert f["mean_intensity"][1] == 10.0
    assert f["elongation"][0] > 5 > f["elongation"][1]
    # a single voxel, or a one-voxel-thick sheet, still has finite elongation
    flat = np.zeros((5, 20, 20), dtype=np.int32)
    flat[2, 5:9, 5:9] = 1
    flat[0, 0, 0] = 2
    g = ins.instance_features(flat, [1, 2], np.ones(flat.shape), ZOOMS)
    assert np.all(np.isfinite(g["elongation"]))
    assert 1.0 <= g["elongation"][1] < 2.0


def test_blockwise_matches_whole(ins):
    """Seed ownership: tiling into overlapping blocks reproduces the whole-volume
    instances exactly, each counted once."""
    rng = np.random.default_rng(3)
    shape = (32, 96, 96)
    ext = np.array(shape) * ZOOMS
    centres = rng.uniform(10, ext - 10, (40, 3))
    img = _puncta(shape, centres, sigma_um=4.0)
    mask = _mask(img)
    args = dict(zooms_um=ZOOMS, sigma_um=3.0, min_dist_um=4.0, h=7.0)

    whole = ins.block_instances(img, mask, (0, 0, 0), shape, (0, 0, 0), **args)

    depth = (6, 12, 12)
    chunk = (16, 48, 48)
    pad = [(d, d) for d in depth]
    img_p = np.pad(img, pad)
    mask_p = np.pad(mask, pad)
    seeds = []
    for z0 in range(0, shape[0], chunk[0]):
        for y0 in range(0, shape[1], chunk[1]):
            for x0 in range(0, shape[2], chunk[2]):
                start = np.array([z0, y0, x0])
                size = np.minimum(chunk, np.array(shape) - start)
                sl = tuple(
                    slice(s, s + n + 2 * d) for s, n, d in zip(start, size, depth)
                )
                res = ins.block_instances(
                    img_p[sl],
                    mask_p[sl],
                    depth,
                    tuple(d + n for d, n in zip(depth, size)),
                    tuple(start - np.array(depth)),
                    **args,
                )
                if res is not None:
                    seeds.append(res["seed"])
    seeds = np.concatenate(seeds)
    assert len(seeds) == len(whole["seed"])
    a = {tuple(s) for s in whole["seed"].astype(int)}
    b = {tuple(s) for s in seeds.astype(int)}
    assert a == b


def test_min_volume_filter(ins):
    shape = (24, 64, 64)
    img = _puncta(shape, [np.array(shape) * ZOOMS / 2])
    mask = _mask(img)
    kw = dict(zooms_um=ZOOMS, sigma_um=3.0, min_dist_um=4.0, h=7.0)
    assert ins.block_instances(img, mask, (0, 0, 0), shape, (0, 0, 0), **kw)
    assert (
        ins.block_instances(
            img, mask, (0, 0, 0), shape, (0, 0, 0), min_volume_um3=1e9, **kw
        )
        is None
    )


def test_h_from_noise(ins):
    # mostly quiet blocks (sd ~2) plus a tail of blocks dense with signal
    sds = [2.0, 2.1, 1.9, 2.2, 2.0, 6.5, 12.0, 17.0]
    h, sd, source = ins.h_from_noise(sds, 3, fallback=7)
    assert source.startswith("noise")
    assert 1.9 <= sd <= 2.2
    assert h == pytest.approx(3 * sd)
    # NaN (too little tissue outside the mask) and zero blocks are ignored
    h, _, source = ins.h_from_noise([np.nan, 0.0, 2.0], 3, fallback=7)
    assert (h, source.startswith("fallback")) == (7.0, True)


def test_block_noise_sd_ignores_mask_and_outside(ins):
    rng = np.random.default_rng(1)
    img = rng.normal(90, 5, (30, 80, 80)).astype(np.float32)
    img[:, :, 60:] = 0  # outside the brain
    mask = np.zeros(img.shape, dtype=bool)
    mask[:, :20, :] = True
    img[mask] += 1000  # signal must not inflate the noise estimate
    sd = ins.block_noise_sd(img, mask, 3.0, 12.0, ZOOMS, 50_000, seed=0)
    # smoothing at 3 um shrinks white noise of sd 5 by roughly a factor of 4-5
    assert 0.3 < sd < 2.5
    # a block that is all mask has no noise estimate
    assert np.isnan(
        ins.block_noise_sd(img, np.ones_like(mask), 3.0, 12.0, ZOOMS, 50_000, 0)
    )


def test_to_physical_matches_zarrnii(ins):
    from zarrnii.analysis import _transform_coordinate_to_physical

    affine = np.array(
        [[0, 0, 0.0035, 1.0], [0, -0.0035, 0, 2.0], [0.005, 0, 0, -3.0], [0, 0, 0, 1]]
    )
    vox = np.array([[1.5, 2.0, 3.0], [10, 20, 30]])
    np.testing.assert_allclose(
        ins.to_physical(vox, affine), _transform_coordinate_to_physical(vox, affine)
    )
