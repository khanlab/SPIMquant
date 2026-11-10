"""Multiscale Frangi fibre orientation analysis (fibre_orientation.py)."""

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

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
def fo():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        spec = spec_from_file_location(
            "fibre_orientation", SCRIPTS_DIR / "fibre_orientation.py"
        )
        module = module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(SCRIPTS_DIR))
    return module


# ----------------------------------------------------------------------
# Phantoms
# ----------------------------------------------------------------------


def unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def fibre_phantom(
    shape, spacing, directions, radius_um=2.0, points=None, background=10.0, seed=0
):
    """Bright Gaussian-profile lines on a constant tissue background.

    ``directions`` are in array-axis (physical) order, one per fibre.  Lines
    pass through ``points`` (physical, array order), random by default.
    """
    rng = np.random.default_rng(seed)
    directions = unit(np.atleast_2d(directions))
    ext = np.asarray(shape) * np.asarray(spacing)
    if points is None:
        points = rng.uniform(0.2, 0.8, size=(len(directions), 3)) * ext
    grids = np.meshgrid(
        *[np.arange(n) * s for n, s in zip(shape, spacing)], indexing="ij"
    )
    coords = np.stack(grids, axis=-1)
    img = np.zeros(shape, dtype=np.float32)
    for d, p in zip(directions, np.atleast_2d(points)):
        r = coords - p
        perp = r - (r @ d)[..., None] * d
        dist2 = np.sum(perp**2, axis=-1)
        img = np.maximum(img, 100.0 * np.exp(-dist2 / (2 * radius_um**2)))
    return (img + background).astype(np.float32)


def angle_deg(a, b):
    """Axial angle between unit vectors (sign-invariant), in degrees."""
    c = np.abs(np.sum(unit(a) * unit(b), axis=-1))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def random_directions(n, rng):
    return unit(rng.normal(size=(n, 3)))


def perturbed_directions(mean, n, spread_deg, rng):
    """Directions scattered around ``mean`` by ~``spread_deg``."""
    mean = unit(mean)
    v = mean + np.tan(np.radians(spread_deg)) * rng.normal(size=(n, 3))
    return unit(v)


SCALES = (1.5, 2.5)
GAMMA = 5.0


# ----------------------------------------------------------------------
# Spherical harmonics
# ----------------------------------------------------------------------


def test_sh_orthonormal(fo):
    # quadrature on a fine Fibonacci sphere
    n = 40000
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    phi = np.pi * (1 + 5**0.5) * i
    r = np.sqrt(1 - z**2)
    v = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)
    y = fo.real_sh_basis(v, 8)
    gram = y.T @ y * (4 * np.pi / n)
    np.testing.assert_allclose(gram, np.eye(y.shape[1]), atol=2e-3)


def test_sh_mrtrix_convention(fo):
    """l=2 terms in closed form (MRtrix3 / dipy tournier07, non-legacy)."""
    rng = np.random.default_rng(0)
    v = random_directions(100, rng)
    x, y, z = v.T
    c = np.sqrt(15 / (4 * np.pi))
    expected = {
        (2, -2): c * x * y,
        (2, -1): -c * y * z,
        (2, 0): np.sqrt(5 / (16 * np.pi)) * (3 * z**2 - 1),
        (2, 1): -c * x * z,
        (2, 2): c / 2 * (x**2 - y**2),
    }
    sh = fo.real_sh_basis(v, 2)
    for (l, m), val in expected.items():
        np.testing.assert_allclose(sh[:, fo.sh_index(l, m)], val, atol=1e-12)


def test_sh_matches_dipy(fo):
    shm = pytest.importorskip("dipy.reconst.shm")
    rng = np.random.default_rng(1)
    v = random_directions(500, rng)
    ref, _, _ = shm.real_sh_tournier(
        8, np.arccos(v[:, 2]), np.arctan2(v[:, 1], v[:, 0]), legacy=False
    )
    np.testing.assert_allclose(fo.real_sh_basis(v, 8), ref, atol=1e-10)


def test_sh_antipodal(fo):
    v = random_directions(50, np.random.default_rng(2))
    np.testing.assert_allclose(fo.real_sh_basis(v, 6), fo.real_sh_basis(-v, 6))


# ----------------------------------------------------------------------
# Eigen-decomposition
# ----------------------------------------------------------------------


def _h6(mats):
    return np.stack(
        [mats[:, i, j] for i, j in [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]], 1
    )


def test_analytic_eigen_matches_lapack(fo):
    rng = np.random.default_rng(10)
    m = rng.normal(size=(20000, 3, 3))
    m = m + m.transpose(0, 2, 1)
    # tube-like: one small, two large negative eigenvalues, random rotation
    q, _ = np.linalg.qr(rng.normal(size=(5000, 3, 3)))
    lam = np.stack(
        [rng.normal(0, 0.01, 5000), -1 - rng.random(5000), -1 - rng.random(5000)], 1
    )
    tubes = np.einsum("nij,nj,nkj->nik", q, lam, q)
    for mats in (m, tubes):
        h6 = _h6(mats).astype(np.float32)
        w_a, v_a = fo.sorted_eigen(h6)
        w_l, v_l = fo.sorted_eigen(h6, method="lapack")
        np.testing.assert_allclose(w_a, w_l, atol=1e-4 * np.abs(w_l).max())
        assert np.percentile(angle_deg(v_a, v_l), 99.9) < 0.1


def test_analytic_eigen_degenerate(fo):
    mats = np.stack(
        [
            np.eye(3),
            np.zeros((3, 3)),
            np.diag([0.0, -1.0, -1.0]),
            np.diag([-1.0, -1.0, 0.5]),
        ]
    )
    w, v = fo.sorted_eigen(_h6(mats).astype(np.float32))
    assert np.all(np.isfinite(w)) and np.all(np.isfinite(v))
    np.testing.assert_allclose(np.linalg.norm(v, axis=1), 1, atol=1e-5)
    np.testing.assert_allclose(np.abs(v[2]), [1, 0, 0], atol=1e-5)
    np.testing.assert_allclose(np.abs(v[3]), [0, 0, 1], atol=1e-5)


# ----------------------------------------------------------------------
# Orientation recovery
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "direction",
    [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 1, 1), (0.3, -1, 0.5)],
)
@pytest.mark.parametrize("spacing", [(1.0, 1.0, 1.0), (2.0, 1.0, 1.0)])
def test_single_fibre_orientation(fo, direction, spacing):
    shape = tuple(int(48 / s) for s in spacing)
    centre = np.asarray(shape) * np.asarray(spacing) / 2
    img = fibre_phantom(shape, spacing, [direction], points=[centre])
    frangi, vec, sidx, tissue = fo.multiscale_frangi(img, spacing, SCALES, GAMMA)
    # voxels within 2 um of the fibre axis, away from the faces
    grids = np.meshgrid(
        *[np.arange(n) * s for n, s in zip(shape, spacing)], indexing="ij"
    )
    r = np.stack(grids, axis=-1) - centre
    d = unit(direction)
    dist = np.linalg.norm(r - (r @ d)[..., None] * d, axis=-1)
    core = np.zeros(shape, bool)
    core[tuple(slice(n // 4, 3 * n // 4) for n in shape)] = True
    on = core & (dist < 2.0)
    assert np.all(frangi[on] > 0)
    err = angle_deg(vec[on], direction)
    assert np.median(err) < 1.5
    assert np.percentile(err, 95) < 3.0


def test_bright_only(fo):
    """A dark line on a bright background gives no response."""
    spacing = (1.0, 1.0, 1.0)
    img = fibre_phantom((32, 32, 32), spacing, [(0, 0, 1)], points=[(16, 16, 16)])
    dark = (img.max() + 10 - img).astype(np.float32)
    frangi, *_ = fo.multiscale_frangi(dark, spacing, SCALES, GAMMA)
    core = (slice(8, 24),) * 3
    assert frangi[core].max() < 1e-3


def test_scale_selection(fo):
    """Thicker fibres select larger scales."""
    spacing = (1.0, 1.0, 1.0)
    scales = (1.0, 2.0, 4.0)
    means = []
    for r in (1.0, 4.0):
        img = fibre_phantom((40, 40, 40), spacing, [(0, 0, 1)], r, [(20, 20, 20)])
        frangi, _, sidx, _ = fo.multiscale_frangi(img, spacing, scales, GAMMA)
        on = frangi > 0.5 * frangi.max()
        means.append(np.mean(np.asarray(scales)[sidx[on]]))
    assert means[1] > means[0]


def test_mask_edge_suppressed(fo):
    """A flat tissue block (no fibres) produces no fibre voxels at its edge."""
    img = np.zeros((40, 40, 40), np.float32)
    img[8:32, 8:32, 8:32] = 50.0
    frangi, *_ = fo.multiscale_frangi(img, (1, 1, 1), SCALES, GAMMA)
    assert frangi.max() < 1e-3


def test_world_rotation(fo):
    # ZYX array with x flipped: array axis 2 -> world -x, axis 0 -> world z
    affine = np.array(
        [[0, 0, -0.5, 1], [0, 2.0, 0, 2], [4.0, 0, 0, 3], [0, 0, 0, 1]], float
    )
    rot = fo.world_rotation(affine)
    np.testing.assert_allclose(rot @ [1, 0, 0], [0, 0, 1])
    np.testing.assert_allclose(rot @ [0, 1, 0], [0, 1, 0])
    np.testing.assert_allclose(rot @ [0, 0, 1], [-1, 0, 0])


def test_downsampled_affine(fo):
    affine = np.diag([2.0, 3.0, 4.0, 1.0])
    affine[:3, 3] = [10, 20, 30]
    out = fo.downsampled_affine(affine, (4, 4, 2))
    np.testing.assert_allclose(np.diag(out)[:3], [8, 12, 8])
    # centre of block 0 is the centre of voxels 0..3 (i.e. index 1.5)
    np.testing.assert_allclose(out[:3, 3], [10 + 2 * 1.5, 20 + 3 * 1.5, 30 + 4 * 0.5])


# ----------------------------------------------------------------------
# Metrics
# ----------------------------------------------------------------------


def _moments_from_vectors(fo, vecs, lmax=4):
    """Moments of one super-voxel holding fibre voxels with ``vecs``."""
    n = len(vecs)
    m = np.zeros(fo.n_moments(lmax))
    m[fo.N_TISSUE] = m[fo.N_FIBRE] = n
    for k, (i, j) in enumerate(fo._TENSOR_IDX):
        m[fo.TENSOR.start + k] = np.sum(vecs[:, i] * vecs[:, j])
    m[fo.SH_START :] = fo.real_sh_basis(vecs, lmax).sum(axis=0)
    return m


def test_metrics_aligned_vs_isotropic(fo):
    rng = np.random.default_rng(3)
    aligned = _moments_from_vectors(fo, np.tile([[0.0, 0.0, 1.0]], (100, 1)))
    iso = _moments_from_vectors(fo, random_directions(20000, rng))
    a = fo.moments_to_metrics(aligned)
    i = fo.moments_to_metrics(iso)
    assert a["odi"] == pytest.approx(0, abs=1e-6)
    assert a["fa"] == pytest.approx(1, abs=1e-6)
    assert abs(a["dirz"]) == pytest.approx(1, abs=1e-6)
    assert i["odi"] == pytest.approx(0.5, abs=0.02)
    assert i["fa"] == pytest.approx(0, abs=0.03)
    assert a["gfa"] > 0.8 and i["gfa"] < 0.1


def test_metrics_planar(fo):
    """Fibres spread in a plane: primary dispersion high, secondary low."""
    rng = np.random.default_rng(4)
    t = rng.uniform(0, np.pi, 5000)
    planar = np.stack([np.cos(t), np.sin(t), np.zeros_like(t)], axis=1)
    m = fo.moments_to_metrics(_moments_from_vectors(fo, planar))
    assert m["odipri"] > 0.45
    assert m["odisec"] < 0.05
    assert abs(m["dirz"]) < 1e-6


def test_dispersion_monotonic(fo):
    rng = np.random.default_rng(5)
    odi = [
        float(
            fo.moments_to_metrics(
                _moments_from_vectors(
                    fo, perturbed_directions([1, 0, 1], 5000, spread, rng)
                )
            )["odi"]
        )
        for spread in (2, 10, 25, 60)
    ]
    assert np.all(np.diff(odi) > 0)


def test_empty_moments_are_zero(fo):
    m = fo.moments_to_metrics(np.zeros((fo.n_moments(6), 2, 2, 2)))
    for name, arr in m.items():
        assert np.all(arr == 0), name
        assert np.all(np.isfinite(arr)), name


def test_aggregate_moments(fo):
    rng = np.random.default_rng(6)
    m = rng.uniform(0, 5, size=(fo.n_moments(2), 4, 4, 4))
    labels = rng.integers(0, 3, size=(4, 4, 4))
    agg = fo.aggregate_moments(m, labels)
    for lab in range(3):
        np.testing.assert_allclose(agg[:, lab], m[:, labels == lab].sum(axis=1))
    assert fo.lmax_from_n_moments(agg.shape[0]) == 2


# ----------------------------------------------------------------------
# Blocks, dask and parameter estimation
# ----------------------------------------------------------------------


def test_reduce_partial_blocks(fo):
    shape = (10, 7, 9)
    rng = np.random.default_rng(7)
    tissue = np.ones(shape, bool)
    fibre = rng.random(shape) > 0.5
    vec = random_directions(np.prod(shape), rng).reshape(shape + (3,))
    m = fo.reduce_to_moments(
        np.ones(shape, np.float32), vec, np.ones(shape), tissue, fibre, (4, 4, 4), 2
    )
    assert m.shape == (fo.n_moments(2), 3, 2, 3)
    assert m[fo.N_TISSUE].sum() == np.prod(shape)
    assert m[fo.N_FIBRE].sum() == fibre.sum()
    assert m[fo.N_TISSUE, 2, 1, 2] == 2 * 3 * 1
    # trace of the tensor sum equals the fibre count (unit vectors)
    trace = m[fo.TENSOR.start] + m[fo.TENSOR.start + 3] + m[fo.TENSOR.start + 5]
    np.testing.assert_allclose(trace, m[fo.N_FIBRE], rtol=1e-5)


def _crossing_phantom(shape=(64, 48, 56), spacing=(1.5, 1.0, 1.0), seed=8):
    rng = np.random.default_rng(seed)
    dirs = np.concatenate(
        [
            perturbed_directions([0, 0, 1], 12, 5, rng),
            perturbed_directions([0, 1, 0], 12, 5, rng),
        ]
    )
    img = fibre_phantom(shape, spacing, dirs, radius_um=1.5, seed=seed)
    # a tissue mask that does not fill the volume
    img[:, :, :6] = 0
    return img, spacing


def test_dask_matches_numpy(fo):
    """Chunked moments equal moments of the whole (zero-padded) volume."""
    da = pytest.importorskip("dask.array")
    img, spacing = _crossing_phantom()
    factor = (4, 4, 4)
    kw = dict(
        spacing_um=spacing,
        scales_um=SCALES,
        gammas=GAMMA,
        threshold=0.2,
        factor=factor,
        lmax=4,
    )
    depth = fo.overlap_depth(SCALES, spacing)
    ref = fo.fibre_moments_block(np.pad(img, [(d, d) for d in depth]), trim=depth, **kw)
    darr = da.from_array(img, chunks=(24, 16, 20))
    with __import__("dask").config.set(scheduler="synchronous"):
        out = fo.fibre_moments(darr, chunk_size=(24, 16, 20), **kw).compute()
    assert out.shape == ref.shape
    assert ref[fo.N_FIBRE].sum() > 100
    np.testing.assert_allclose(out, ref, rtol=1e-4, atol=1e-3)


def test_crossing_fibres_odf(fo):
    """Two populations: dispersion is planar and the ODF has two peaks."""
    img, spacing = _crossing_phantom()
    depth = fo.overlap_depth(SCALES, spacing)
    m = fo.fibre_moments_block(
        np.pad(img, [(d, d) for d in depth]),
        spacing,
        SCALES,
        GAMMA,
        0.2,
        img.shape,  # one super-voxel for the whole volume
        lmax=8,
        trim=depth,
    )
    met = fo.moments_to_metrics(m[:, 0, 0, 0])
    # primary dispersion (in the y-z plane) much larger than secondary
    assert met["odipri"] > 0.35
    assert met["odisec"] < 0.15
    # identity rotation: fibres along array axes 2 and 1 are world z and y,
    # so the ODF peaks there and is small along world x (array axis 0)
    coeffs = fo.odf_coefficients(m[:, 0, 0, 0])
    amp = fo.real_sh_basis(np.eye(3), 8) @ coeffs
    assert amp[2] > 3 * amp[0] and amp[1] > 3 * amp[0]


def test_rotation_applied(fo):
    """Fibres along array axis 0 map to the world axis given by the affine."""
    spacing = (1.0, 1.0, 1.0)
    img = fibre_phantom((32, 32, 32), spacing, [(1, 0, 0)], points=[(16, 16, 16)])
    rot = fo.world_rotation(
        np.array([[0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1.0]])
    )
    depth = fo.overlap_depth(SCALES, spacing)
    m = fo.fibre_moments_block(
        np.pad(img, [(d, d) for d in depth]),
        spacing,
        SCALES,
        GAMMA,
        0.2,
        (32, 32, 32),
        lmax=2,
        rotation=rot,
        trim=depth,
    )
    met = fo.moments_to_metrics(m[:, 0, 0, 0])
    assert abs(met["dirz"]) > 0.99


def test_parameter_estimation(fo):
    da = pytest.importorskip("dask.array")
    img, spacing = _crossing_phantom()
    darr = da.from_array(img, chunks=(32, 24, 28))
    blocks = fo.tissue_blocks(darr, 4)
    assert len(blocks) == 4
    est = fo.sample_frangi_parameters(
        darr, spacing, SCALES, blocks, n_samples_per_block=5000
    )
    assert est["gammas"].shape == (len(SCALES),)
    assert np.all(est["gammas"] > 0)
    thr = est["threshold"]
    # the threshold separates on-fibre voxels from the background
    frangi, *_ = fo.multiscale_frangi(img, spacing, SCALES, est["gammas"])
    near = img > 60
    far = (img > 0) & (img < 10.5)
    assert np.mean(frangi[near] > thr) > 0.7
    assert np.mean(frangi[far] > thr) < 0.05


def test_metrics_from_dask(fo):
    da = pytest.importorskip("dask.array")
    rng = np.random.default_rng(9)
    m = np.abs(rng.normal(size=(fo.n_moments(2), 6, 6, 6))) + 1
    out = fo.metrics_from_dask_moments(da.from_array(m, chunks=(4, 3, 3, 3)))
    np.testing.assert_allclose(out.compute(), fo.metrics_array(m), rtol=1e-6)


@pytest.mark.parametrize("orientation", ["RAS", "LPI"])
def test_zarrnii_wrapper(fo, orientation):
    zarrnii = pytest.importorskip("zarrnii")
    da = pytest.importorskip("dask.array")
    img, spacing = _crossing_phantom()
    znimg = zarrnii.ZarrNii.from_darr(
        da.from_array(img[None], chunks=(1, 32, 24, 28)),
        spacing=tuple(s / 1000 for s in spacing),  # mm
        origin=(1.0, 2.0, 3.0),
        orientation=orientation,
    )
    out, params = fo.zarrnii_fibre_moments(
        znimg, SCALES, 4, lmax=4, chunk_size=(32, 24, 28), n_sample_blocks=4
    )
    assert out.darr.shape == (fo.n_moments(4), 16, 12, 14)
    np.testing.assert_allclose(
        out.affine.matrix,
        fo.downsampled_affine(znimg.affine.matrix, (4, 4, 4)),
        atol=1e-9,
    )
    m = out.darr.compute()
    assert m[fo.N_FIBRE].sum() > 100
    assert params["threshold"] > 0


def test_store_moments_by_slab(fo, tmp_path):
    da = pytest.importorskip("dask.array")
    img, spacing = _crossing_phantom()
    darr = da.from_array(img, chunks=(16, 24, 28))
    kw = dict(
        spacing_um=spacing,
        scales_um=SCALES,
        gammas=GAMMA,
        threshold=0.2,
        factor=(4, 4, 4),
        lmax=2,
        chunk_size=(16, 24, 28),
    )
    lazy = fo.fibre_moments(darr, **kw)
    assert len(lazy.chunks[1]) == 4
    stored = fo.store_moments_by_slab(
        fo.fibre_moments(darr, **kw), str(tmp_path / "m.zarr"), slab_chunks=3
    )
    with __import__("dask").config.set(scheduler="synchronous"):
        np.testing.assert_allclose(stored.compute(), lazy.compute(), rtol=1e-6)


# ----------------------------------------------------------------------
# PSF equalisation and soft threshold
# ----------------------------------------------------------------------


def _blobs(shape, spacing, sigma_um, n=200, seed=11):
    """Gaussian blobs of per-axis sigma ``sigma_um`` (array order)."""
    rng = np.random.default_rng(seed)
    coords = np.stack(
        np.meshgrid(*[np.arange(k) * s for k, s in zip(shape, spacing)], indexing="ij"),
        -1,
    )
    ext = np.asarray(shape) * np.asarray(spacing)
    img = np.zeros(shape, np.float32)
    for c in rng.uniform(0.1, 0.9, size=(n, 3)) * ext:
        img += 100 * np.exp(-np.sum(((coords - c) / sigma_um) ** 2, -1) / 2)
    return img + 10 + rng.normal(0, 0.5, shape).astype(np.float32)


def test_psf_equalizing_sigma(fo):
    fwhm = np.array([6.0, 3.6, 3.6])
    extra = fo.psf_equalizing_sigma(fwhm)
    assert extra[0] == 0
    total = (fwhm / 2.3548) ** 2 + extra**2
    np.testing.assert_allclose(total, total[0])


def test_psf_correction_removes_axial_bias(fo):
    """Point-like blobs blurred by an axially elongated PSF look like axial
    tubes; equalising the PSF makes them blob-like again, so their Frangi
    response collapses (and they fall below a fibre threshold)."""
    spacing = (2.0, 1.625, 1.625)
    psf_fwhm = np.array([6.0, 3.6, 3.6])
    img = _blobs((96, 128, 128), spacing, psf_fwhm / 2.3548, n=60)
    core = img > 60
    res = []
    for extra in (None, fo.psf_equalizing_sigma(psf_fwhm)):
        fr, vec, *_ = fo.multiscale_frangi(
            img, spacing, (2.0, 3.0), 50.0, psf_sigma_um=extra
        )
        sel = core & (fr > 0)
        res.append((np.abs(vec[sel][:, 0]).mean(), fr[core].mean()))
    assert res[0][0] > 0.85  # spurious axial orientation without correction
    assert res[1][1] < res[0][1] / 3  # blobs no longer tube-like


def test_puncta_fwhm(fo):
    spacing = (2.0, 1.625, 1.625)
    fwhm = np.array([6.0, 3.5, 3.5])
    img = _blobs((80, 96, 96), spacing, fwhm / 2.3548, n=60)
    rows = fo.puncta_fwhm(img, spacing, quantile=0.99)
    assert len(rows) > 20
    np.testing.assert_allclose(np.median(rows, 0), fwhm, rtol=0.15)


def test_soft_weights(fo):
    f = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.8])
    np.testing.assert_allclose(fo.soft_weights(f, 0.4, 1.0), [0, 0, 0, 0, 0, 1])
    np.testing.assert_allclose(
        fo.soft_weights(f, 0.4, 0.25), [0, 0, 1 / 3, 2 / 3, 1, 1], atol=1e-6
    )


def test_soft_dask_matches_numpy(fo):
    da = pytest.importorskip("dask.array")
    img, spacing = _crossing_phantom()
    kw = dict(
        spacing_um=spacing,
        scales_um=SCALES,
        gammas=GAMMA,
        threshold=0.4,
        factor=(4, 4, 4),
        lmax=2,
        soft_lo=0.25,
        psf_sigma_um=(1.0, 0.0, 0.0),
    )
    depth = fo.overlap_depth(SCALES, spacing, psf_sigma_um=kw["psf_sigma_um"])
    ref = fo.fibre_moments_block(np.pad(img, [(d, d) for d in depth]), trim=depth, **kw)
    with __import__("dask").config.set(scheduler="synchronous"):
        out = fo.fibre_moments(
            da.from_array(img, chunks=(24, 16, 20)), chunk_size=(24, 16, 20), **kw
        ).compute()
    np.testing.assert_allclose(out, ref, rtol=1e-4, atol=1e-3)
    # soft weights: fractional fibre counts, more than with the hard threshold
    hard = fo.fibre_moments_block(
        np.pad(img, [(d, d) for d in depth]), trim=depth, **{**kw, "soft_lo": 1.0}
    )
    assert ref[fo.N_FIBRE].sum() > hard[fo.N_FIBRE].sum()
