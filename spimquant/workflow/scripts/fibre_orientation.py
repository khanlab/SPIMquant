"""Multiscale Frangi fibre orientation analysis (FOA) for large 3D images.

Estimates local fibre orientations with a multiscale Hessian (Frangi)
filter and summarises them over super-voxel blocks as orientation
distribution functions (ODFs, real spherical harmonics) and orientation
tensors, from which dispersion and anisotropy metrics are derived.

The method follows Foa3D (Sorelli et al. 2023, Sci Rep 13:4160;
https://github.com/lens-biophotonics/Foa3D, MIT licence), adapted for
chunked (dask) processing:

* Gaussian scales are given in micrometres and derivatives are taken in
  physical units, so anisotropic voxels give correct angles without
  resampling to isotropic.
* The Frangi background parameter ``gamma`` and the fibre-mask threshold
  are global (estimated once from sampled blocks), not per slice/chunk,
  so results do not depend on the chunking.
* The image is blurred with a mask-normalised Gaussian (zeros are treated
  as missing, not as dark tissue), so the brain-mask edge does not look
  like a fibre.
* Per-voxel orientations are never written out: each chunk is reduced
  directly to additive per-block *moments* (voxel counts, orientation
  tensor sums and spherical-harmonic sums).  Because they are sums, blocks
  can be pooled exactly into coarser blocks or atlas regions before the
  metrics are computed.

Layout (kept free of snakemake so it can move into zarrnii later):

1. numpy kernels on in-memory blocks (Hessian, Frangi, SH, reductions),
2. metrics from moments (numpy, any leading shape),
3. dask drivers (``fibre_moments``, ``sample_frangi_parameters``),
4. a thin ZarrNii wrapper (``zarrnii_fibre_moments``).

Array conventions: spatial arrays are in array-axis order (e.g. ZYX for
zarrnii) with ``spacing_um`` in the same order.  Orientation vectors are
computed in that frame and rotated to world (RAS, x/y/z) coordinates with
the 3x3 direction-cosine matrix of the image affine (``world_rotation``)
before the ODF is formed, so the SH coefficients are in the MRtrix3
convention and can be viewed directly in ``mrview``.
"""

import math

import numpy as np
from scipy import ndimage as ndi
from scipy.special import gammaln, lpmv

# Gaussian kernels are truncated at this many sigmas; overlap depth follows.
TRUNCATE = 4.0

# Voxels per batch in the analytic eigensolver (bounds its float64 temporaries).
_EIGEN_BATCH = 1 << 18

# Below this, the blur weight G(M) is treated as "no tissue nearby".
_MIN_WEIGHT = 1e-3

# Moment layout (axis 0 of moment arrays); SH sums follow at SH_START.
N_TISSUE = 0  # tissue voxels in the block
N_FIBRE = 1  # fibre voxels in the block
SUM_FRANGI = 2  # sum of the Frangi response over tissue voxels
SUM_SCALE = 3  # sum of the selected scale (um) over fibre voxels
TENSOR = slice(4, 10)  # sum of v v^T over fibre voxels: xx xy xz yy yz zz
SH_START = 10
_TENSOR_IDX = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))

METRIC_NAMES = (
    "fibredensity",  # fibre voxels / tissue voxels
    "frangi",  # mean Frangi response over tissue
    "scale",  # mean selected Frangi scale over fibre voxels (um)
    "odi",  # total orientation dispersion index (0 aligned .. 1)
    "odipri",  # dispersion in the primary plane
    "odisec",  # dispersion in the secondary plane
    "odianis",  # anisotropy of the dispersion
    "fa",  # fractional anisotropy of the orientation tensor
    "gfa",  # generalised fractional anisotropy of the ODF
    "dirx",  # principal fibre direction (world x, y, z; sign arbitrary)
    "diry",
    "dirz",
)


# ----------------------------------------------------------------------
# Parameters and geometry
# ----------------------------------------------------------------------


def sigma_vox(sigma_um, spacing_um):
    """Per-axis Gaussian sigma in voxels for a sigma in micrometres."""
    return tuple(float(sigma_um) / float(s) for s in spacing_um)


def overlap_depth(scales_um, spacing_um, truncate=TRUNCATE):
    """Per-axis overlap (voxels) a block needs for the largest scale.

    The kernel radius plus 2 voxels for the two central differences.
    """
    smax = max(float(s) for s in scales_um)
    return tuple(
        int(math.ceil(truncate * sv)) + 2 for sv in sigma_vox(smax, spacing_um)
    )


def n_sh_coeffs(lmax):
    """Number of even-order real SH coefficients up to ``lmax``."""
    if lmax % 2 or lmax < 0:
        raise ValueError(f"lmax must be even and >= 0, got {lmax}")
    return (lmax + 1) * (lmax + 2) // 2


def n_moments(lmax):
    """Length of the moment axis for a given ``lmax``."""
    return SH_START + n_sh_coeffs(lmax)


def lmax_from_n_moments(n):
    """Inverse of :func:`n_moments`."""
    ncoef = n - SH_START
    lmax = int(round((-3 + math.sqrt(1 + 8 * ncoef)) / 2))
    if n_sh_coeffs(lmax) != ncoef:
        raise ValueError(f"{n} is not a valid number of moments")
    return lmax


def world_rotation(affine):
    """Direction cosines (3x3) of an index-to-world affine.

    Column ``j`` is the world (x, y, z) unit vector of array axis ``j``,
    so a vector ``v`` with components along the array axes (in physical
    units) is ``world_rotation(affine) @ v`` in world coordinates.  Shear
    is not supported (the columns are only normalised).
    """
    a = np.asarray(affine, dtype=np.float64)[:3, :3]
    return a / np.linalg.norm(a, axis=0, keepdims=True)


def downsampled_affine(affine, factor):
    """Affine of a grid whose voxels are ``factor`` blocks of the input.

    Voxel centres of the block grid sit at the centres of the blocks.
    """
    a = np.array(affine, dtype=np.float64)
    f = np.asarray(factor, dtype=np.float64)
    out = a.copy()
    out[:3, :3] = a[:3, :3] * f[None, :]
    out[:3, 3] = a[:3, 3] + a[:3, :3] @ ((f - 1) / 2)
    return out


# ----------------------------------------------------------------------
# Hessian and Frangi (numpy, one block)
# ----------------------------------------------------------------------


def normalized_blur(img, weight, sigma, truncate=TRUNCATE):
    """Gaussian blur of ``img`` using only voxels where ``weight`` is set.

    Returns ``(blur, den)`` where ``den = G(weight)``; ``blur`` is zero
    where ``den`` is negligible.
    """
    num = ndi.gaussian_filter(img * weight, sigma, truncate=truncate)
    den = ndi.gaussian_filter(weight, sigma, truncate=truncate)
    out = np.zeros_like(num)
    np.divide(num, den, out=out, where=den > _MIN_WEIGHT)
    return out, den


def _iter_hessian_components(smoothed, spacing_um):
    """Yield the six Hessian components one at a time (see below)."""
    grads = {}
    for i, j in _TENSOR_IDX:
        if i == j:
            kernel = np.array([1.0, -2.0, 1.0], dtype=np.float32) / spacing_um[i] ** 2
            yield ndi.correlate1d(smoothed, kernel, axis=i, mode="nearest")
        else:
            if i not in grads:
                grads[i] = np.gradient(smoothed, spacing_um[i], axis=i)
            yield np.gradient(grads[i], spacing_um[j], axis=j)
            if j == 2:  # last mixed term using grads[i]
                del grads[i]


def hessian_components(smoothed, spacing_um):
    """Six Hessian components (00, 01, 02, 11, 12, 22) in physical units.

    Diagonal terms use the compact second difference ``[1, -2, 1] / h^2``
    and mixed terms nested central differences.  (Foa3D nests central
    differences for all terms, i.e. a ``[1, 0, -2, 0, 1] / 4h^2`` stencil
    on the diagonal, which over-smooths coarsely sampled axes and biases
    angles by several degrees on anisotropic voxels.)
    """
    return list(_iter_hessian_components(smoothed, spacing_um))


def _eigvals_sym3(a, b, c, d, e, f):
    """Closed-form eigenvalues of symmetric 3x3 matrices (ascending).

    ``[[a, b, c], [b, d, e], [c, e, f]]`` given element-wise (float64), by
    the trigonometric method (Smith 1961).
    """
    q = (a + d + f) / 3
    p1 = b * b + c * c + e * e
    a0, d0, f0 = a - q, d - q, f - q
    p = np.sqrt((a0 * a0 + d0 * d0 + f0 * f0 + 2 * p1) / 6)
    safe = np.where(p > 0, p, 1.0)
    det = (
        a0 * (d0 * f0 - e * e) - b * (b * f0 - e * c) + c * (b * e - d0 * c)
    ) / safe**3
    phi = np.arccos(np.clip(det / 2, -1.0, 1.0)) / 3
    hi = q + 2 * p * np.cos(phi)
    lo = q + 2 * p * np.cos(phi + 2 * np.pi / 3)
    return np.stack([lo, 3 * q - hi - lo, hi], axis=1)


def sorted_eigen(h6, method="analytic"):
    """Eigen-decomposition of symmetric 3x3 matrices given as (N, 6).

    Returns ``(eigvals, vec)``: eigenvalues (N, 3) sorted by increasing
    absolute value, and the unit eigenvector (N, 3) of the smallest-
    magnitude eigenvalue, i.e. the fibre axis.

    ``method="analytic"`` (default) uses closed-form eigenvalues and a
    cross product for the eigenvector, ~10x faster than LAPACK; matrices
    where that eigenvector is ill-defined (repeated smallest eigenvalue)
    fall back to ``np.linalg.eigh``, as does ``method="lapack"``.
    """
    h6 = np.asarray(h6)
    n = h6.shape[0]
    if method == "lapack":
        return _sorted_eigen_lapack(h6)
    if n > _EIGEN_BATCH:
        # bounded float64 temporaries (~150 B per voxel per batch)
        w = np.empty((n, 3), dtype=np.float32)
        vec = np.empty((n, 3), dtype=np.float32)
        for start in range(0, n, _EIGEN_BATCH):
            sl = slice(start, start + _EIGEN_BATCH)
            w[sl], vec[sl] = sorted_eigen(h6[sl], method)
        return w, vec
    h = h6.astype(np.float64)
    a, b, c, d, e, f = (h[:, k] for k in range(6))
    w = _eigvals_sym3(a, b, c, d, e, f)
    order = np.argsort(np.abs(w), axis=1)
    w = np.take_along_axis(w, order, axis=1)
    lam = w[:, 0]
    r0 = np.stack([a - lam, b, c], axis=1)
    r1 = np.stack([b, d - lam, e], axis=1)
    r2 = np.stack([c, e, f - lam], axis=1)
    cands = np.stack([np.cross(r0, r1), np.cross(r0, r2), np.cross(r1, r2)])
    norms = np.linalg.norm(cands, axis=2)
    best = np.argmax(norms, axis=0)
    vec = cands[best, np.arange(n)]
    nrm = norms[best, np.arange(n)]
    # ill-conditioned: the null space of H - lam I is not 1D
    scale = np.max(np.abs(w), axis=1)
    bad = nrm <= 1e-6 * np.maximum(scale, 1e-30) ** 2
    vec = vec / np.where(bad, 1.0, nrm)[:, None]
    if bad.any():
        wb, vb = _sorted_eigen_lapack(h6[bad])
        w[bad] = wb
        vec[bad] = vb
    return w.astype(np.float32), vec.astype(np.float32)


def _sorted_eigen_lapack(h6):
    n = h6.shape[0]
    mat = np.empty((n, 3, 3), dtype=np.float64)
    for k, (i, j) in enumerate(_TENSOR_IDX):
        mat[:, i, j] = h6[:, k]
        mat[:, j, i] = h6[:, k]
    w, v = np.linalg.eigh(mat)
    order = np.argsort(np.abs(w), axis=1)
    w = np.take_along_axis(w, order, axis=1)
    vec = np.take_along_axis(v, order[:, None, :1], axis=2)[..., 0]
    return w.astype(np.float32), vec.astype(np.float32)


def frangi_vesselness(eigvals, alpha, beta, gamma):
    """Frangi tubularness for bright tubes from sorted eigenvalues (..., 3).

    ``eigvals`` are sorted by increasing absolute value.  Voxels that are
    not bright tubes (l2 > 0 or l3 > 0) score zero.
    """
    l1, l2, l3 = (np.abs(eigvals[..., k]) for k in range(3))
    with np.errstate(divide="ignore", invalid="ignore"):
        ra = np.where(l3 > 0, l2 / l3, 0.0)
        l23 = np.sqrt(l2 * l3)
        rb = np.where(l23 > 0, l1 / l23, 0.0)
    s2 = l1**2 + l2**2 + l3**2
    v = (
        (1 - np.exp(-(ra**2) / (2 * alpha**2)))
        * np.exp(-(rb**2) / (2 * beta**2))
        * (1 - np.exp(-s2 / (2 * float(gamma) ** 2)))
    )
    bright = (eigvals[..., 1] < 0) & (eigvals[..., 2] < 0)
    return np.where(bright & np.isfinite(v), v, 0.0).astype(np.float32)


def structureness(eigvals):
    """Frobenius norm of the Hessian, ``sqrt(l1^2 + l2^2 + l3^2)``."""
    return np.sqrt(np.sum(np.square(eigvals, dtype=np.float64), axis=-1))


def _scale_eigen(
    img, weight, candidates, spacing_um, scale_um, edge_min_weight, tubes_only
):
    """Scale-normalised Hessian eigen at one scale.

    Analyses voxels in ``candidates`` that have at least ``edge_min_weight``
    tissue in their Gaussian neighbourhood and, if ``tubes_only``, a
    negative Hessian trace (necessary for a bright tube: l2, l3 < 0 and
    |l1| <= |l2|).  Returns ``(idx, eigvals, vec)`` with ``idx`` the flat
    indices of the analysed voxels.
    """
    blur, den = normalized_blur(img, weight, sigma_vox(scale_um, spacing_um))
    idx = np.flatnonzero(candidates & (den >= edge_min_weight))
    del den
    if idx.size == 0:
        return idx, np.zeros((0, 3), np.float32), np.zeros((0, 3), np.float32)
    # gather each component at the analysed voxels as soon as it is computed,
    # rather than holding all six full blocks
    h6 = np.empty((idx.size, 6), dtype=np.float32)
    for k, comp in enumerate(_iter_hessian_components(blur, spacing_um)):
        h6[:, k] = comp.ravel()[idx]
    del blur
    h6 *= np.float32(float(scale_um) ** 2)
    if tubes_only:
        keep = (h6[:, 0] + h6[:, 3] + h6[:, 5]) < 0
        idx, h6 = idx[keep], h6[keep]
    w, vec = sorted_eigen(h6)
    return idx, w, vec


def _core_mask(shape, trim):
    core = np.zeros(shape, dtype=bool)
    core[tuple(slice(t, n - t) for t, n in zip(trim, shape))] = True
    return core


def multiscale_frangi(
    img,
    spacing_um,
    scales_um,
    gammas,
    alpha=0.001,
    beta=1.0,
    edge_min_weight=0.5,
    trim=None,
):
    """Multiscale Frangi filter of one 3D block.

    Parameters
    ----------
    img : ndarray (3D)
        Intensity; voxels <= 0 are treated as outside the tissue.
    spacing_um : sequence of 3 float
        Voxel size per array axis, in micrometres.
    scales_um : sequence of float
        Gaussian scales (sigma) in micrometres.
    gammas : float or sequence of float
        Frangi background parameter, one per scale (or shared).
    alpha, beta : float
        Frangi plate and blob parameters (Foa3D defaults: 0.001, 1).
    edge_min_weight : float
        A voxel is only analysed at a scale if at least this fraction of
        its Gaussian neighbourhood is tissue (suppresses the mask edge).
    trim : sequence of 3 int, optional
        Overlap to drop from each face: only the core is analysed and
        returned.

    Returns
    -------
    frangi : ndarray float32
        Maximum response over scales (0 outside tissue).
    vec : ndarray float32 (... , 3)
        Unit fibre axis at the best scale, components along the array
        axes (physical frame); zero where ``frangi`` is 0.
    scale_idx : ndarray int8
        Index of the best scale, -1 where ``frangi`` is 0.
    tissue : ndarray bool
    """
    img = np.asarray(img, dtype=np.float32)
    tissue = img > 0
    weight = tissue.astype(np.float32)
    candidates = tissue
    if trim is not None and any(trim):
        candidates = tissue & _core_mask(img.shape, trim)
    gammas = np.broadcast_to(np.asarray(gammas, dtype=np.float64), (len(scales_um),))

    best = np.zeros(img.size, dtype=np.float32)
    vec = np.zeros((img.size, 3), dtype=np.float32)
    scale_idx = np.full(img.size, -1, dtype=np.int8)
    for s, scale_um in enumerate(scales_um):
        idx, w, v = _scale_eigen(
            img, weight, candidates, spacing_um, scale_um, edge_min_weight, True
        )
        if idx.size == 0:
            continue
        resp = frangi_vesselness(w, alpha, beta, gammas[s])
        better = resp > best[idx]
        sel = idx[better]
        best[sel] = resp[better]
        vec[sel] = v[better]
        scale_idx[sel] = s
    out = (
        best.reshape(img.shape),
        vec.reshape(img.shape + (3,)),
        scale_idx.reshape(img.shape),
        tissue,
    )
    if trim is not None and any(trim):
        core = tuple(slice(t, n - t) for t, n in zip(trim, img.shape))
        out = tuple(o[core] for o in out)
    return out


def sample_eigenvalues(
    img, spacing_um, scales_um, n_samples, edge_min_weight=0.5, seed=0, trim=None
):
    """Sorted Hessian eigenvalues at random tissue voxels of a block.

    The same voxels are used for every scale, so the multiscale Frangi
    response can be evaluated afterwards for any ``gamma``.

    Parameters
    ----------
    trim : sequence of 3 int, optional
        Only sample voxels at least this far from the block faces (the
        overlap added around a chunk).

    Returns
    -------
    ndarray (n_scales, n, 3)
        NaN where a voxel is too close to the mask edge at that scale.
    """
    img = np.asarray(img, dtype=np.float32)
    tissue = img > 0
    weight = tissue.astype(np.float32)
    candidates = tissue
    if trim is not None and any(trim):
        candidates = tissue & _core_mask(img.shape, trim)
    cand = np.flatnonzero(candidates)
    rng = np.random.default_rng(seed)
    if cand.size > n_samples:
        cand = np.sort(rng.choice(cand, size=n_samples, replace=False))
    out = np.full((len(scales_um), cand.size, 3), np.nan, dtype=np.float32)
    if cand.size == 0:
        return out
    sampled = np.zeros(img.size, dtype=bool)
    sampled[cand] = True
    sampled = sampled.reshape(img.shape)
    for s, scale_um in enumerate(scales_um):
        idx, w, _ = _scale_eigen(
            img, weight, sampled, spacing_um, scale_um, edge_min_weight, False
        )
        # idx is a sorted subset of cand (voxels too near the mask edge drop)
        out[s, np.searchsorted(cand, idx)] = w
    return out


def estimate_gammas(eig_samples, quantile=0.999, factor=0.5):
    """Global Frangi ``gamma`` per scale from sampled eigenvalues.

    Foa3D uses half the maximum structureness of each slice; here a high
    quantile over the whole sample is used instead, which is global (no
    chunk seams) and robust to isolated outliers.
    """
    gammas = []
    for w in eig_samples:
        s = structureness(w)
        s = s[np.isfinite(s)]
        g = factor * float(np.quantile(s, quantile)) if s.size else 0.0
        gammas.append(g if g > 0 else 1.0)
    return np.asarray(gammas)


def sampled_frangi(eig_samples, gammas, alpha=0.001, beta=1.0):
    """Multiscale (max) Frangi response at sampled voxels."""
    resp = [
        frangi_vesselness(np.nan_to_num(w, nan=0.0), alpha, beta, g)
        for w, g in zip(eig_samples, gammas)
    ]
    return np.max(np.stack(resp), axis=0)


def estimate_threshold(frangi_values, method="yen"):
    """Global fibre-mask threshold on the Frangi response.

    ``method`` is ``"yen"`` (as Foa3D), ``"otsu"`` or ``"li"``, applied to
    the *positive* sampled responses (most tissue voxels are not tube-like
    and score exactly 0; that spike would otherwise pull the threshold to
    ~0), or a number, used as the threshold itself.
    """
    try:
        return float(method)
    except (TypeError, ValueError):
        pass
    from skimage import filters

    values = np.asarray(frangi_values, dtype=np.float64)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size < 2 or values.max() <= values.min():
        return float("inf")
    fn = {
        "yen": filters.threshold_yen,
        "otsu": filters.threshold_otsu,
        "li": filters.threshold_li,
    }[method]
    return float(fn(values))


# ----------------------------------------------------------------------
# Spherical harmonics (MRtrix3 convention)
# ----------------------------------------------------------------------


def sh_index(l, m):
    """Column of the even-order real SH coefficient (l, m)."""
    return l * (l + 1) // 2 + m


def real_sh_basis(vecs_xyz, lmax):
    """Even-order real SH basis evaluated at unit vectors.

    Uses the orthonormal MRtrix3 convention (with the Condon-Shortley
    phase, as in scipy):
    ``Y_lm = sqrt(2) Im(Y_l^|m|)`` for m < 0, ``Y_l0`` for m = 0 and
    ``sqrt(2) Re(Y_l^m)`` for m > 0, ordered by l then m.  Matches dipy's
    ``tournier07`` basis with ``legacy=False``.

    Parameters
    ----------
    vecs_xyz : ndarray (N, 3)
        Directions in world (x, y, z); need not be normalised.

    Returns
    -------
    ndarray float64 (N, n_sh_coeffs(lmax))
    """
    v = np.asarray(vecs_xyz, dtype=np.float64).reshape(-1, 3)
    norm = np.linalg.norm(v, axis=1)
    norm[norm == 0] = 1.0
    cos_t = np.clip(v[:, 2] / norm, -1.0, 1.0)
    phi = np.arctan2(v[:, 1], v[:, 0])
    out = np.empty((v.shape[0], n_sh_coeffs(lmax)), dtype=np.float64)
    for l in range(0, lmax + 1, 2):
        for m in range(0, l + 1):
            nlm = math.sqrt(
                (2 * l + 1)
                / (4 * math.pi)
                * math.exp(gammaln(l - m + 1) - gammaln(l + m + 1))
            )
            # scipy's lpmv includes the Condon-Shortley phase, as MRtrix3 does
            p = nlm * lpmv(m, l, cos_t)
            if m == 0:
                out[:, sh_index(l, 0)] = p
            else:
                out[:, sh_index(l, m)] = math.sqrt(2) * p * np.cos(m * phi)
                out[:, sh_index(l, -m)] = math.sqrt(2) * p * np.sin(m * phi)
    return out


# ----------------------------------------------------------------------
# Block reduction to moments
# ----------------------------------------------------------------------


def _block_ids(shape, factor):
    """Flat output-block index of every voxel, and the output grid shape."""
    grid = tuple(-(-n // f) for n, f in zip(shape, factor))
    axes = [np.arange(n) // f for n, f in zip(shape, factor)]
    ids = np.ravel_multi_index(np.meshgrid(*axes, indexing="ij"), grid)
    return ids.ravel(), grid


def reduce_to_moments(frangi, vec_world, scale_um, tissue, fibre, factor, lmax):
    """Reduce per-voxel results of a block to per-super-voxel moments.

    Parameters
    ----------
    frangi : ndarray (3D)
        Frangi response.
    vec_world : ndarray (3D + (3,))
        Unit fibre axes in world (x, y, z).
    scale_um : ndarray (3D)
        Selected scale (um) per voxel.
    tissue, fibre : ndarray bool (3D)
        Tissue mask and fibre mask (fibre implies tissue).
    factor : sequence of 3 int
        Super-voxel size in voxels per axis.  Blocks at the far faces may
        be partial.
    lmax : int
        Maximum even SH order.

    Returns
    -------
    ndarray float32 (n_moments(lmax), *grid)
    """
    ids, grid = _block_ids(frangi.shape, factor)
    nb = int(np.prod(grid))
    out = np.zeros((n_moments(lmax), nb), dtype=np.float64)
    t = tissue.ravel()
    f = fibre.ravel()
    out[N_TISSUE] = np.bincount(ids, weights=t, minlength=nb)
    out[N_FIBRE] = np.bincount(ids, weights=f, minlength=nb)
    out[SUM_FRANGI] = np.bincount(ids[t], weights=frangi.ravel()[t], minlength=nb)
    fid = ids[f]
    if fid.size:
        out[SUM_SCALE] = np.bincount(fid, weights=scale_um.ravel()[f], minlength=nb)
        v = vec_world.reshape(-1, 3)[f].astype(np.float64)
        for k, (i, j) in enumerate(_TENSOR_IDX):
            out[TENSOR.start + k] = np.bincount(
                fid, weights=v[:, i] * v[:, j], minlength=nb
            )
        sh = real_sh_basis(v, lmax)
        for k in range(sh.shape[1]):
            out[SH_START + k] = np.bincount(fid, weights=sh[:, k], minlength=nb)
    return out.reshape((-1,) + grid).astype(np.float32)


def fibre_moments_block(
    img,
    spacing_um,
    scales_um,
    gammas,
    threshold,
    factor,
    lmax=6,
    rotation=None,
    alpha=0.001,
    beta=1.0,
    edge_min_weight=0.5,
    trim=None,
):
    """Moments of one (possibly overlapped) block.

    ``trim`` is the overlap to remove from each face after filtering; the
    trimmed shape is what ``factor`` divides into super-voxels.
    """
    frangi, vec, sidx, tissue = multiscale_frangi(
        img, spacing_um, scales_um, gammas, alpha, beta, edge_min_weight, trim
    )
    fibre = tissue & (frangi > threshold) & (sidx >= 0)
    rot = np.eye(3) if rotation is None else np.asarray(rotation, np.float64)
    vec_world = vec @ rot.T.astype(np.float32)
    scale_lut = np.append(np.asarray(scales_um, np.float32), np.float32(0))
    scale_um = scale_lut[sidx]  # -1 -> 0
    return reduce_to_moments(frangi, vec_world, scale_um, tissue, fibre, factor, lmax)


# ----------------------------------------------------------------------
# Metrics from moments (numpy; moment axis first, any trailing shape)
# ----------------------------------------------------------------------


def _safe_div(a, b):
    out = np.zeros(np.broadcast(a, b).shape, dtype=np.float64)
    np.divide(a, b, out=out, where=b > 0)
    return out


def orientation_tensor(moments):
    """Mean orientation tensor (..., 3, 3) (trace 1 where there are fibres)."""
    m = np.asarray(moments, dtype=np.float64)
    n = m[N_FIBRE]
    t = np.zeros(n.shape + (3, 3))
    for k, (i, j) in enumerate(_TENSOR_IDX):
        t[..., i, j] = t[..., j, i] = _safe_div(m[TENSOR.start + k], n)
    return t


def dispersion_indices(evals_desc):
    """Foa3D dispersion indices from tensor eigenvalues (..., 3), descending.

    With l1 >= l2 >= l3, each index is ``1 - (2/pi) arctan(l1 / x)``:
    total (x = sqrt(l2 l3)), primary (x = l2), secondary (x = l3) and
    anisotropy (x = |l2 - l3|).  0 means perfectly aligned fibres.
    """
    l1, l2, l3 = (np.abs(evals_desc[..., k]) for k in range(3))
    c = 2.0 / np.pi
    return {
        "odi": 1 - c * np.arctan2(l1, np.sqrt(l2 * l3)),
        "odipri": 1 - c * np.arctan2(l1, l2),
        "odisec": 1 - c * np.arctan2(l1, l3),
        "odianis": 1 - c * np.arctan2(l1, np.abs(l2 - l3)),
    }


def tensor_fa(evals):
    """Fractional anisotropy from eigenvalues (..., 3)."""
    md = evals.mean(axis=-1, keepdims=True)
    num = np.sqrt(np.sum((evals - md) ** 2, axis=-1))
    den = np.sqrt(np.sum(evals**2, axis=-1))
    return np.sqrt(1.5) * _safe_div(num, den)


def odf_coefficients(moments):
    """Mean SH coefficients (n_coeffs, ...) of the fibre-orientation ODF."""
    m = np.asarray(moments, dtype=np.float64)
    return _safe_div(m[SH_START:], m[N_FIBRE][None])


def sh_gfa(coeffs):
    """Generalised fractional anisotropy from SH coefficients (axis 0)."""
    c = np.asarray(coeffs, dtype=np.float64)
    power = np.sum(c**2, axis=0)
    return np.sqrt(np.clip(1 - _safe_div(c[0] ** 2, power), 0, 1)) * (power > 0)


def moments_to_metrics(moments):
    """All scalar/vector metrics from moments.

    Parameters
    ----------
    moments : ndarray (n_moments, ...)

    Returns
    -------
    dict name -> ndarray (...), keys as in ``METRIC_NAMES``.  Voxels or
    regions without fibres get 0 for every orientation metric.
    """
    m = np.asarray(moments, dtype=np.float64)
    has = m[N_FIBRE] > 0
    t = orientation_tensor(m)
    w, v = np.linalg.eigh(t)
    w = w[..., ::-1]  # descending
    disp = dispersion_indices(w)
    out = {
        "fibredensity": _safe_div(m[N_FIBRE], m[N_TISSUE]),
        "frangi": _safe_div(m[SUM_FRANGI], m[N_TISSUE]),
        "scale": _safe_div(m[SUM_SCALE], m[N_FIBRE]),
        **{k: np.where(has, d, 0.0) for k, d in disp.items()},
        "fa": np.where(has, tensor_fa(w), 0.0),
        "gfa": sh_gfa(odf_coefficients(m)),
    }
    e1 = v[..., :, 2]
    for k, name in enumerate(("dirx", "diry", "dirz")):
        out[name] = np.where(has, e1[..., k], 0.0)
    return {k: out[k].astype(np.float32) for k in METRIC_NAMES}


def metrics_array(moments):
    """``moments_to_metrics`` stacked along axis 0 (``METRIC_NAMES`` order)."""
    d = moments_to_metrics(moments)
    return np.stack([d[k] for k in METRIC_NAMES])


def aggregate_moments(moments, labels, n_labels=None):
    """Sum moments over label regions.

    Parameters
    ----------
    moments : ndarray (n_moments, *grid)
    labels : ndarray int (*grid)
        Region label per super-voxel (0 = background is kept as a row).

    Returns
    -------
    ndarray float64 (n_moments, n_labels); feed to ``moments_to_metrics``.
    """
    lab = np.asarray(labels).ravel().astype(np.int64)
    if n_labels is None:
        n_labels = int(lab.max()) + 1 if lab.size else 1
    m = np.asarray(moments, dtype=np.float64).reshape(moments.shape[0], -1)
    return np.stack([np.bincount(lab, weights=row, minlength=n_labels) for row in m])


# ----------------------------------------------------------------------
# Dask drivers
# ----------------------------------------------------------------------


def aligned_chunks(shape, chunk_size, factor):
    """Chunks along each axis that are multiples of ``factor``.

    Every chunk except the last has a size that is the multiple of the
    factor nearest to the requested ``chunk_size`` (at least one factor).
    """
    out = []
    for n, c, f in zip(shape, chunk_size, factor):
        c = max(int(f), int(round(c / f)) * int(f))
        full, rem = divmod(int(n), c)
        out.append((c,) * full + ((rem,) if rem else ()))
    return tuple(out)


def _moments_chunk(block, depth=None, **kwargs):
    return fibre_moments_block(block, trim=depth, **kwargs)


def fibre_moments(
    darr,
    spacing_um,
    scales_um,
    gammas,
    threshold,
    factor,
    lmax=6,
    rotation=None,
    alpha=0.001,
    beta=1.0,
    edge_min_weight=0.5,
    chunk_size=None,
):
    """Per-super-voxel fibre moments of a 3D dask array.

    The array is rechunked so chunk sizes are multiples of ``factor`` (so
    no super-voxel straddles chunks), overlapped by the largest kernel
    radius (zero padding at the array boundary, which counts as
    non-tissue), filtered and reduced chunk by chunk.

    Returns
    -------
    dask.array (n_moments(lmax), ceil(shape / factor)) float32
    """
    import dask.array as da

    factor = tuple(int(f) for f in factor)
    if chunk_size is None:
        chunk_size = tuple(c[0] for c in darr.chunks)
    depth = overlap_depth(scales_um, spacing_um)
    chunks = aligned_chunks(darr.shape, chunk_size, factor)
    # overlap needs every chunk to be at least as deep as the overlap
    for ax, (ch, d) in enumerate(zip(chunks, depth)):
        if min(ch) < d and len(ch) > 1:
            raise ValueError(
                f"axis {ax}: chunk {min(ch)} smaller than overlap {d}; "
                "increase chunk_size"
            )
    darr = darr.rechunk(chunks)
    ov = da.overlap.overlap(
        darr, depth=dict(enumerate(depth)), boundary={0: 0, 1: 0, 2: 0}
    )
    out_chunks = ((n_moments(lmax),),) + tuple(
        tuple(-(-c // f) for c in ch) for ch, f in zip(chunks, factor)
    )
    return da.map_blocks(
        _moments_chunk,
        ov,
        depth=depth,
        spacing_um=tuple(float(s) for s in spacing_um),
        scales_um=tuple(float(s) for s in scales_um),
        gammas=np.broadcast_to(np.asarray(gammas, float), (len(scales_um),)),
        threshold=float(threshold),
        factor=factor,
        lmax=lmax,
        rotation=None if rotation is None else np.asarray(rotation, float),
        alpha=alpha,
        beta=beta,
        edge_min_weight=edge_min_weight,
        new_axis=0,
        chunks=out_chunks,
        dtype=np.float32,
        meta=np.empty((0, 0, 0, 0), dtype=np.float32),
    )


def metrics_from_dask_moments(moments):
    """``metrics_array`` applied chunkwise to a dask moment array."""
    import dask.array as da

    moments = moments.rechunk({0: -1})
    return da.map_blocks(
        metrics_array,
        moments,
        chunks=((len(METRIC_NAMES),),) + moments.chunks[1:],
        dtype=np.float32,
    )


def sample_frangi_parameters(
    darr,
    spacing_um,
    scales_um,
    block_ids,
    n_samples_per_block=20000,
    alpha=0.001,
    beta=1.0,
    edge_min_weight=0.5,
    gamma_quantile=0.999,
    gamma_factor=0.5,
    threshold_method="yen",
    seed=0,
):
    """Estimate global ``gammas`` and fibre threshold from sampled chunks.

    Parameters
    ----------
    darr : dask.array (3D)
    block_ids : list of tuple
        Chunk indices (of ``darr.chunks``) to sample, e.g. tissue-rich
        chunks chosen from a low-resolution level.

    Returns
    -------
    dict with ``gammas`` (per scale), ``threshold``, ``frangi_samples``
    and ``eig_samples`` (for QC).
    """
    import dask
    import dask.array as da

    depth = overlap_depth(scales_um, spacing_um)
    # same zero padding as fibre_moments; blocks are passed to the tasks as
    # dask arrays (no compute inside a task)
    ov = da.overlap.overlap(
        darr, depth=dict(enumerate(depth)), boundary={0: 0, 1: 0, 2: 0}
    )
    parts = dask.compute(
        *[
            dask.delayed(sample_eigenvalues)(
                ov.blocks[tuple(b)],
                spacing_um,
                scales_um,
                n_samples_per_block,
                edge_min_weight=edge_min_weight,
                seed=seed + k,
                trim=depth,
            )
            for k, b in enumerate(block_ids)
        ]
    )
    parts = [p for p in parts if p.shape[1]]
    if not parts:
        raise ValueError("No tissue voxels in the sampled blocks")
    eig = np.concatenate(parts, axis=1)
    gammas = estimate_gammas(eig, quantile=gamma_quantile, factor=gamma_factor)
    frangi = sampled_frangi(eig, gammas, alpha, beta)
    threshold = estimate_threshold(frangi, threshold_method)
    return dict(
        gammas=gammas, threshold=threshold, frangi_samples=frangi, eig_samples=eig
    )


def tissue_blocks(darr, n_blocks, min_fraction=0.5, seed=0):
    """Pick up to ``n_blocks`` random chunks with at least ``min_fraction``
    nonzero voxels (falls back to the most tissue-rich chunks).

    Reads the whole array once; for very large images pass ``block_ids``
    chosen from a low-resolution level instead.
    """
    import dask.array as da

    frac = da.map_blocks(
        lambda b: np.array((b > 0).mean(), ndmin=b.ndim),
        darr,
        chunks=tuple((1,) * len(c) for c in darr.chunks),
        dtype=np.float64,
    ).compute()
    ids = list(zip(*np.nonzero(frac > 0)))
    rng = np.random.default_rng(seed)
    rich = [b for b in ids if frac[b] >= min_fraction]
    rng.shuffle(rich)
    chosen = rich[:n_blocks]
    if len(chosen) < n_blocks:
        rest = sorted(
            (b for b in ids if frac[b] < min_fraction), key=lambda b: -frac[b]
        )
        chosen += rest[: n_blocks - len(chosen)]
    return [tuple(int(i) for i in b) for b in chosen]


# ----------------------------------------------------------------------
# ZarrNii wrapper
# ----------------------------------------------------------------------


def zarrnii_fibre_moments(
    znimg,
    scales_um,
    factor,
    lmax=6,
    gammas=None,
    threshold=None,
    block_ids=None,
    n_sample_blocks=16,
    chunk_size=(128, 128, 128),
    spacing_scale=1000.0,
    **kwargs,
):
    """Fibre moments of a single-channel ZarrNii, as a ZarrNii.

    Parameters
    ----------
    znimg : ZarrNii
        Single-channel image (channel axis first, spatial axes after).
    scales_um : sequence of float
        Frangi scales (sigma) in micrometres.
    factor : int or sequence of 3 int
        Super-voxel size in voxels (per spatial axis).
    gammas, threshold : optional
        Global Frangi parameters; estimated from sampled chunks if None.
    block_ids : list of tuple, optional
        Chunks to sample for the estimate (default: random tissue chunks).
    spacing_scale : float
        Multiplier from ``znimg.get_zooms()`` units to micrometres (1000
        for zooms in mm, as in spimquant).
    **kwargs
        Passed to ``fibre_moments`` / ``sample_frangi_parameters``
        (``alpha``, ``beta``, ``edge_min_weight``, ...).

    Returns
    -------
    (ZarrNii, dict)
        Moments with one channel per moment on the super-voxel grid, and
        the parameters used (``gammas``, ``threshold``, ...).
    """
    darr = znimg.darr
    if darr.ndim != 4 or darr.shape[0] != 1:
        raise ValueError(f"Expected a single-channel image, got {darr.shape}")
    darr = darr[0]
    factor = (int(factor),) * 3 if np.isscalar(factor) else tuple(map(int, factor))
    spacing_um = tuple(float(z) * spacing_scale for z in znimg.get_zooms())
    affine = np.asarray(
        znimg.affine.matrix if hasattr(znimg.affine, "matrix") else znimg.affine
    )
    rotation = world_rotation(affine)

    filt_keys = ("alpha", "beta", "edge_min_weight")
    filt = {k: kwargs[k] for k in filt_keys if k in kwargs}
    params = {}
    if gammas is None or threshold is None:
        sample_arr = darr.rechunk(aligned_chunks(darr.shape, chunk_size, factor))
        if block_ids is None:
            block_ids = tissue_blocks(sample_arr, n_sample_blocks)
        est_keys = filt_keys + (
            "n_samples_per_block",
            "gamma_quantile",
            "gamma_factor",
            "threshold_method",
            "seed",
        )
        params = sample_frangi_parameters(
            sample_arr,
            spacing_um,
            scales_um,
            block_ids,
            **{k: kwargs[k] for k in est_keys if k in kwargs},
        )
        gammas = params["gammas"] if gammas is None else gammas
        threshold = params["threshold"] if threshold is None else threshold
    params.update(
        gammas=np.broadcast_to(gammas, (len(scales_um),)), threshold=threshold
    )

    moments = fibre_moments(
        darr,
        spacing_um,
        scales_um,
        gammas,
        threshold,
        factor,
        lmax=lmax,
        rotation=rotation,
        chunk_size=chunk_size,
        **filt,
    )
    return _moments_zarrnii(znimg, moments, factor), params


def _moments_zarrnii(znimg, moments, factor):
    """ZarrNii holding ``moments`` on the super-voxel grid of ``znimg``.

    Voxel centres sit at block centres, i.e. the affine is
    ``downsampled_affine(znimg.affine, factor)``; orientation metadata is
    kept.  The half-block shift of the translation follows the sign of each
    axis in the affine (zarrnii 0.22's ``downsample`` always shifts in +,
    which misplaces flipped axes, e.g. LPI, by ``f - 1`` input voxels).
    """
    from zarrnii import ZarrNii

    ngff = znimg.ngff_image
    spatial = ["z", "y", "x"] if znimg.axes_order == "ZYX" else ["x", "y", "z"]
    a = np.asarray(getattr(znimg.affine, "matrix", znimg.affine))[:3, :3]
    signs = np.sign(a[np.argmax(np.abs(a), axis=0), np.arange(3)])
    scale = dict(ngff.scale)
    translation = dict(ngff.translation or {})
    for dim, f, sgn in zip(spatial, factor, signs):
        translation[dim] = translation.get(dim, 0.0) + sgn * scale[dim] * (f - 1) / 2
        scale[dim] *= f
    try:
        from zarrnii.core import _derive_ngff_image

        new = _derive_ngff_image(
            ngff, data=moments, scale=scale, translation=translation
        )
    except ImportError:  # pragma: no cover - older/newer zarrnii
        import ngff_zarr as nz

        new = nz.NgffImage(
            data=moments,
            dims=ngff.dims,
            scale=scale,
            translation=translation,
            name="fibre_moments",
        )
    return ZarrNii(
        ngff_image=new,
        axes_order=znimg.axes_order,
        xyz_orientation=znimg.xyz_orientation,
    )
