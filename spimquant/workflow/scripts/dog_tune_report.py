"""DoG segmentation tuning report.

Helps choose ``s`` and ``k`` of ``dog+s{s}r{r}k{k}`` by translating each
setting into quantities that mean something biologically, measured on a
sample of tissue blocks from the subject:

* **Contrast detection limit.**  Synthetic puncta (3D Gaussians) are added
  at known contrast -- peak height relative to the local background
  intensity -- and the contrast at which 50% and 90% of them are detected is
  found for every (s, k).  The DoG is linear in the image (the brain mask is
  fixed), so the response of every punctum at every contrast follows from
  two DoGs per block: of the image and of unit-contrast puncta.
* **Noise false-discovery rate.**  The same blocks are rebuilt as a smooth
  background plus white noise at the measured, intensity-dependent noise
  level, and segmented with the same thresholds; detections there are false
  positives due to noise.
* **Object density** (connected components per mm^3 of tissue) and the
  masked volume fraction for every (s, k).
* **Scale census.**  Sizes of the strongest blobs in the image, from maxima
  of the scale-normalised Laplacian of Gaussian over position and scale.
* **Montages** of the masks for a few (s, k) on three blocks.

k units (the robust sd of the DoG) are computed exactly as dogthresh does:
same chunking, same blocks and samples, so the k values in this report give
the same threshold in a ``dogthresh`` run of the same image.
"""

import base64
import io
import json

import numpy as np
from dogthresh import (
    _TRUNCATE,
    ellipsoid,
    normalized_blur,
    robust_threshold,
    sigma_vox,
    upsample_nearest_block,
)
from scipy import ndimage as ndi

_RANDOM_STATE = 7

# A 3D Gaussian blob of sd b gives the largest -s^2 * (in-plane Laplacian)
# response at s = b * sqrt(2/3)
_CENSUS_TO_BLOB_SD = 1.0 / np.sqrt(2.0 / 3.0)

# 6-neighbour Laplacian; the sum of squares of its weights is 42, so for white
# noise of sd n its response has sd n * sqrt(42) (Immerkaer's estimator)
_LAPLACIAN = np.zeros((3, 3, 3), dtype=np.float32)
_LAPLACIAN[1, 1, 1] = -6
for _axis in range(3):
    for _side in (0, 2):
        _idx = [1, 1, 1]
        _idx[_axis] = _side
        _LAPLACIAN[tuple(_idx)] = 1


# --------------------------------------------------------------------------- #
# Per-block analysis
# --------------------------------------------------------------------------- #


def weighted_dog(img, weight, sigma1, sigma2):
    """DoG of ``img`` blurred over ``weight`` only, zero outside ``weight``.

    Equals ``dogthresh.dog_response`` when ``weight`` is ``img > 0``, and is
    linear in ``img`` for a fixed weight.
    """
    g1 = normalized_blur(img, weight, sigma1)
    g2 = normalized_blur(img, weight, sigma2)
    return np.where(weight > 0, g1 - g2, 0).astype(np.float32), g1


def noise_model(img, tissue, zooms_um, bg_sigma_um=12.0, n_bins=10, min_voxels=1000):
    """Voxel noise sd as a function of the local background intensity.

    The noise sd is estimated from the 6-neighbour Laplacian (which removes
    smooth structure) as ``1.4826 * MAD / sqrt(42)``, in bins of the local
    background (``bg_sigma_um`` blur).  Returns ``(bg, centres, sds)``, the
    background image and the per-bin model, or ``(bg, None, None)`` when
    there is too little tissue.
    """
    lap = ndi.convolve(img, _LAPLACIAN, mode="nearest") / np.sqrt(42.0)
    bg = normalized_blur(
        img, tissue.astype(np.float32), sigma_vox(bg_sigma_um, zooms_um)
    )
    sel = ndi.binary_erosion(tissue, iterations=2)
    if sel.sum() < min_voxels * 2:
        return bg, None, None
    b, l = bg[sel], lap[sel]
    edges = np.unique(np.percentile(b, np.linspace(1, 99, n_bins + 1)))
    centres, sds = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (b >= lo) & (b < hi)
        if m.sum() >= min_voxels:
            centres.append(float(np.median(b[m])))
            sds.append(1.4826 * float(np.median(np.abs(l[m] - np.median(l[m])))))
    if not centres:
        return bg, None, None
    return bg, np.asarray(centres), np.asarray(sds)


def noise_null_image(img, tissue, zooms_um, rng):
    """A structureless image with the measured noise.

    Constant intensity (the tissue median) plus white noise whose sd follows
    ``noise_model`` at the local background level of the real image.  A
    constant mean gives the DoG nothing to respond to, so everything detected
    in it is due to noise.  Returns ``(null_image, noise_sd_at_median)``.
    """
    bg, centres, sds = noise_model(img, tissue, zooms_um)
    if centres is None:
        return None, None
    noise_sd = np.interp(bg, centres, sds).astype(np.float32)
    level = float(np.median(img[tissue]))
    noise = rng.standard_normal(img.shape, dtype=np.float32) * noise_sd
    null = np.where(tissue, level + noise, 0).astype(np.float32)
    return null, float(np.interp(level, centres, sds))


def scale_norm_log(img, scale_um, zooms_um):
    """``-s^2`` x the in-plane Laplacian of the 3D Gaussian-smoothed image.

    The z derivative is left out: with z voxels several micrometres thick it
    is dominated by sampling at small scales.
    """
    sv = sigma_vox(scale_um, zooms_um)
    lap = sum(
        ndi.gaussian_filter(img, sv, order=[2 if a == ax else 0 for a in range(3)])
        / zooms_um[ax] ** 2
        for ax in (1, 2)
    )
    return (-(scale_um**2) * lap).astype(np.float32)


def scale_census(img, valid, zooms_um, scales_um, n_top, core=None):
    """Blob sd (um) of the ``n_top`` strongest scale-space blobs in ``valid``.

    Blobs are local maxima of the scale-normalised LoG over position and
    scale; maxima at the smallest or largest scale are dropped, since their
    true scale is outside the range.  With ``core`` (a tuple of slices), the
    LoG is computed on all of ``img`` but only ``img[core]`` is searched, and
    ``valid`` has the shape of the core.
    """
    core = core or tuple(slice(None) for _ in img.shape)
    stack = np.stack([scale_norm_log(img, s, zooms_um)[core] for s in scales_um])
    peak = (stack == ndi.maximum_filter(stack, size=3)) & (stack > 0) & valid[None]
    peak[0] = peak[-1] = False
    idx = np.nonzero(peak)
    vals = stack[idx]
    order = np.argsort(vals)[::-1][:n_top]
    return np.asarray(scales_um)[idx[0][order]] * _CENSUS_TO_BLOB_SD


def place_spikes(allowed, zooms_um, n, min_sep_um, rng, max_tries=20000):
    """Up to ``n`` voxel positions in ``allowed``, at least ``min_sep_um`` apart."""
    cand = np.argwhere(allowed)
    if len(cand) == 0:
        return np.zeros((0, 3), dtype=int)
    zooms = np.asarray(zooms_um)
    pts = []
    for i in rng.permutation(len(cand))[:max_tries]:
        p = cand[i] * zooms
        if not pts or np.min(np.linalg.norm(np.asarray(pts) - p, axis=1)) >= min_sep_um:
            pts.append(p)
            if len(pts) == n:
                break
    return np.round(np.asarray(pts).reshape(-1, 3) / zooms).astype(int)


def spike_image(shape, pts, amps, sd_um, zooms_um):
    """Sum of 3D Gaussians of sd ``sd_um`` and peak ``amps`` at voxel ``pts``."""
    out = np.zeros(shape, dtype=np.float32)
    zooms = np.asarray(zooms_um)
    half = np.ceil(4 * sd_um / zooms).astype(int)
    for p, a in zip(pts, amps):
        lo = np.maximum(p - half, 0)
        hi = np.minimum(p + half + 1, shape)
        g = [
            np.exp(-(((np.arange(l, h) - c) * z) ** 2) / (2 * sd_um**2))
            for l, h, c, z in zip(lo, hi, p, zooms)
        ]
        out[tuple(slice(l, h) for l, h in zip(lo, hi))] += (
            a * g[0][:, None, None] * g[1][None, :, None] * g[2][None, None, :]
        )
    return out


def analyse_block(
    block,
    pool_up,
    margin,
    zooms_um,
    stats,
    ratio,
    k_grid,
    montage_k,
    spike_sds_um,
    n_spikes,
    spike_sep_um,
    pool_thresh,
    pool_dilate_um,
    census_scales_um,
    census_n_top,
    seed,
):
    """All tuning measurements for one overlapped block.

    ``block`` is the corrected image with ``margin`` voxels of context on
    each side, ``pool_up`` the low-resolution pool mask upsampled onto it, and
    ``stats`` maps each s (um) to the dogthresh ``(median, robust_sd)``.
    Measurements are made in the block core only.
    """
    rng = np.random.default_rng(seed)
    img = np.asarray(block, dtype=np.float32)
    tissue = img > 0
    weight = tissue.astype(np.float32)
    core = tuple(slice(m, n - m) for m, n in zip(margin, img.shape))
    pool_struct = ellipsoid(pool_dilate_um, zooms_um)
    pool_up = np.asarray(pool_up, dtype=bool)
    eroded = ndi.binary_erosion(tissue, iterations=3)
    vox_vol_mm3 = float(np.prod(zooms_um)) * 1e-9

    out = dict(s={}, tissue_fraction=float(tissue[core].mean()))

    null, noise_sd = noise_null_image(img, tissue, zooms_um, rng)
    out["noise_sd"] = noise_sd
    bg12 = normalized_blur(img, weight, sigma_vox(12.0, zooms_um))
    mid = img.shape[0] // 2
    slab = slice(max(mid - 1, 0), mid + 2)
    out["montage_img"] = img[slab].max(axis=0)[core[1:]]
    out["montage"] = {}

    any_pool = np.zeros(img.shape, dtype=bool)
    for s_um, (med, sd) in stats.items():
        s1 = sigma_vox(s_um, zooms_um)
        s2 = sigma_vox(s_um * ratio, zooms_um)
        dog, g1 = weighted_dog(img, weight, s1, s2)

        # pool exclusion exactly as dogthresh.segment_block
        pool = pool_up & (g1 > pool_thresh)
        if pool_struct is not None and pool.any():
            pool = ndi.binary_dilation(pool, structure=pool_struct)
        any_pool |= pool
        valid = tissue & ~pool
        valid_core = valid[core]
        valid_mm3 = float(valid_core.sum()) * vox_vol_mm3

        res = dict(valid_mm3=valid_mm3, n_objects=[], volume_fraction=[], n_null=[])
        dog_core = dog[core]
        for k in k_grid:
            pos = (dog_core > med + k * sd) & valid_core
            res["n_objects"].append(int(ndi.label(pos)[1]))
            res["volume_fraction"].append(float(pos.sum() / max(valid_core.sum(), 1)))
        for k in montage_k:
            pos = (dog > med + k * sd) & valid
            out["montage"][(s_um, k)] = pos[slab].any(axis=0)[core[1:]]

        if null is not None:
            dnull_core = weighted_dog(null, weight, s1, s2)[0][core]
            for k in k_grid:
                pos = (dnull_core > med + k * sd) & valid_core
                res["n_null"].append(int(ndi.label(pos)[1]))
        else:
            res["n_null"] = None

        # puncta on background: in eroded, pool-free tissue of the core and
        # away from anything already above k = 3
        near = ndi.binary_dilation(
            dog > med + 3 * sd, structure=ellipsoid(spike_sep_um / 3, zooms_um)
        )
        allowed = np.zeros(img.shape, dtype=bool)
        allowed[core] = (eroded & valid & ~near)[core]
        res["spikes"] = {}
        for sd_um in spike_sds_um:
            pts = place_spikes(allowed, zooms_um, n_spikes, spike_sep_um, rng)
            if len(pts) == 0:
                res["spikes"][sd_um] = (np.zeros(0), np.zeros(0))
                continue
            amps = bg12[tuple(pts.T)]
            spikes = spike_image(img.shape, pts, amps, sd_um, zooms_um) * weight
            dspk = weighted_dog(spikes, weight, s1, s2)[0]
            # detection score in k units at contrast c: base + c * gain
            base = (dog[tuple(pts.T)] - med) / sd
            gain = dspk[tuple(pts.T)] / sd
            res["spikes"][sd_um] = (base.astype(np.float32), gain.astype(np.float32))
        out["s"][s_um] = res

    out["census_blob_sd"] = scale_census(
        img, (tissue & ~any_pool)[core], zooms_um, census_scales_um, census_n_top, core
    )
    return out


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #


def contrast_at_recall(base, gain, k, recall):
    """Contrast at which a fraction ``recall`` of the puncta is detected at k.

    A punctum is detected at contrast c when ``base + c * gain > k``; puncta
    with no positive gain are never detected.  Returns inf when fewer than
    ``recall`` of the puncta can be detected at any contrast.
    """
    base, gain = np.asarray(base), np.asarray(gain)
    if base.size == 0:
        return float("nan")
    c_min = np.where(gain > 0, (k - base) / np.where(gain > 0, gain, 1), np.inf)
    c_min = np.maximum(c_min, 0)
    return float(np.quantile(c_min, recall, method="higher"))


def recall_at(base, gain, k, contrast):
    """Fraction of puncta detected at ``contrast`` and threshold ``k``."""
    base, gain = np.asarray(base), np.asarray(gain)
    if base.size == 0:
        return float("nan")
    return float(np.mean(base + contrast * gain > k))


def k_for_target(base, gain, contrast, recall):
    """Largest k that still detects a fraction ``recall`` at ``contrast``."""
    base, gain = np.asarray(base), np.asarray(gain)
    if base.size == 0:
        return float("nan")
    return float(np.quantile(base + contrast * gain, 1 - recall, method="lower"))


def method_string(s_um, ratio, k):
    def fmt(v):
        return f"{v:g}".replace(".", "p")

    return f"dog+s{fmt(s_um)}r{fmt(ratio)}k{fmt(k)}"


def summarise(blocks, stats, ratio, k_grid, spike_sds_um, target):
    """Pool per-block results into per-(s, k) curves and suggestions."""
    summary = {}
    for s_um, (med, sd) in stats.items():
        per = [b["s"][s_um] for b in blocks]
        valid_mm3 = sum(p["valid_mm3"] for p in per)
        n_obj = np.sum([p["n_objects"] for p in per], axis=0)
        vf = np.average(
            [p["volume_fraction"] for p in per],
            axis=0,
            weights=[max(p["valid_mm3"], 1e-12) for p in per],
        )
        nulls = [p["n_null"] for p in per if p["n_null"] is not None]
        n_null = np.sum(nulls, axis=0) if nulls else None
        spikes = {}
        for sd_um in spike_sds_um:
            base = np.concatenate([p["spikes"][sd_um][0] for p in per])
            gain = np.concatenate([p["spikes"][sd_um][1] for p in per])
            spikes[sd_um] = dict(
                n=int(base.size),
                c50=[contrast_at_recall(base, gain, k, 0.5) for k in k_grid],
                c90=[contrast_at_recall(base, gain, k, 0.9) for k in k_grid],
                _base=base,
                _gain=gain,
            )
        tb = spikes[target["spike_sd_um"]]
        k_sens = k_for_target(
            tb["_base"], tb["_gain"], target["contrast"], target["recall"]
        )
        # noise floor: the smallest k on the grid whose noise FDR is within
        # the limit (FDR falls with k); the suggestion is never below it
        fdr = (
            [float(nn / max(no, 1)) for nn, no in zip(n_null, n_obj)]
            if n_null is not None
            else None
        )
        k_noise = None
        if fdr is not None:
            ok = [k for k, f in zip(k_grid, fdr) if f <= target["max_noise_fdr"]]
            k_noise = float(ok[0]) if ok else None
        # lowest k at which the noise FDR is known to be within the limit
        k_low = k_noise if k_noise is not None else float(min(k_grid))
        k_sugg = np.floor(k_sens * 2) / 2 if np.isfinite(k_sens) else float("nan")
        if np.isfinite(k_sugg) and k_sugg >= k_low:
            status = "meets target"
        else:
            status = (
                f"target not met: needs k ≤ {k_sens:.1f}, below k = {k_low:g}"
                if np.isfinite(k_sens)
                else "target not met"
            )
            k_sugg = float("nan")
        summary[s_um] = dict(
            median=med,
            robust_sd=sd,
            threshold=[med + k * sd for k in k_grid],
            valid_mm3=valid_mm3,
            objects_per_mm3=(n_obj / max(valid_mm3, 1e-12)).tolist(),
            volume_fraction=vf.tolist(),
            noise_fdr=fdr,
            noise_objects_per_mm3=(
                (n_null / max(valid_mm3, 1e-12)).tolist()
                if n_null is not None
                else None
            ),
            spikes=spikes,
            k_for_target=k_sens,
            k_noise_floor=k_noise,
            k_suggested=float(k_sugg),
            status=status,
            method=method_string(s_um, ratio, k_sugg) if np.isfinite(k_sugg) else None,
        )
        if np.isfinite(k_sugg):
            summary[s_um]["at_suggested"] = {
                str(sd_um): dict(
                    c50=contrast_at_recall(v["_base"], v["_gain"], k_sugg, 0.5),
                    c90=contrast_at_recall(v["_base"], v["_gain"], k_sugg, 0.9),
                )
                for sd_um, v in spikes.items()
            }
    return summary


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #


def _png(fig):
    import matplotlib.pyplot as plt

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _figures(summary, k_grid, spike_sds_um, census, blocks, montage_ids, montage_k):
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    s_list = list(summary)
    cmap = plt.get_cmap("viridis")
    colors = {s: cmap(i / max(len(s_list) - 1, 1)) for i, s in enumerate(s_list)}
    figs = {}

    # detection limits
    fig, axes = plt.subplots(
        1, len(spike_sds_um), figsize=(5.2 * len(spike_sds_um), 4), squeeze=False
    )
    for ax, sd_um in zip(axes[0], spike_sds_um):
        for s in s_list:
            sp = summary[s]["spikes"][sd_um]
            ax.plot(
                k_grid,
                np.array(sp["c90"]) * 100,
                "-o",
                ms=3,
                color=colors[s],
                label=f"s = {s:g} µm",
            )
            ax.plot(k_grid, np.array(sp["c50"]) * 100, "--", color=colors[s], lw=1)
        ax.set_yscale("log")
        ax.set_ylim(10, 1000)
        ax.set_xlabel("k (robust sd)")
        ax.set_ylabel("contrast above local background (%)")
        ax.set_title(f"puncta of sd {sd_um:g} µm: solid 90%, dashed 50% detected")
        ax.grid(alpha=0.3, which="both")
    axes[0][0].legend(fontsize=8)
    figs["limits"] = _png(fig)

    # density and noise FDR
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for s in s_list:
        sm = summary[s]
        axes[0].plot(
            k_grid,
            sm["objects_per_mm3"],
            "-o",
            ms=3,
            color=colors[s],
            label=f"s = {s:g} µm",
        )
        axes[1].plot(
            k_grid, np.array(sm["volume_fraction"]) * 100, "-o", ms=3, color=colors[s]
        )
        if sm["noise_fdr"] is not None:
            axes[2].plot(
                k_grid, np.maximum(sm["noise_fdr"], 1e-4), "-o", ms=3, color=colors[s]
            )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("objects per mm³ of tissue")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("masked volume (% of tissue)")
    axes[2].set_yscale("log")
    axes[2].set_ylabel("noise FDR (floored at 1e-4)")
    fig.subplots_adjust(wspace=0.35)
    for ax in axes:
        ax.set_xlabel("k (robust sd)")
        ax.grid(alpha=0.3, which="both")
    axes[0].legend(fontsize=8)
    figs["counts"] = _png(fig)

    # census
    fig, ax = plt.subplots(figsize=(6, 3.5))
    if census.size:
        vals, counts = np.unique(np.round(census, 2), return_counts=True)
        ax.bar(range(len(vals)), counts, color="steelblue")
        ax.set_xticks(range(len(vals)), [f"{v:.1f}" for v in vals])
        med = float(np.median(census))
        ax.axvline(
            float(np.interp(med, vals, range(len(vals)))),
            color="k",
            lw=1,
            label=f"median {med:.1f} µm",
        )
        ax.legend(fontsize=8)
    ax.set_xlabel("blob sd (µm)")
    ax.set_ylabel("blobs")
    ax.set_title("sizes of the strongest blobs")
    figs["census"] = _png(fig)

    # montages: rows s, columns k
    figs["montage"] = []
    for bi in montage_ids:
        b = blocks[bi]
        im = b["montage_img"]
        lo, hi = np.percentile(im[im > 0], [1, 99.5]) if (im > 0).any() else (0, 1)
        fig, axes = plt.subplots(
            len(s_list),
            len(montage_k) + 1,
            figsize=(2.3 * (len(montage_k) + 1), 2.3 * len(s_list)),
            squeeze=False,
        )
        for i, s in enumerate(s_list):
            axes[i][0].imshow(im, cmap="gray", vmin=lo, vmax=hi)
            axes[i][0].set_ylabel(f"s = {s:g} µm")
            for j, k in enumerate(montage_k):
                ax = axes[i][j + 1]
                ax.imshow(im, cmap="gray", vmin=lo, vmax=hi)
                m = np.ma.masked_where(~b["montage"][(s, k)], b["montage"][(s, k)])
                ax.imshow(m, cmap="autumn", alpha=0.55, interpolation="nearest")
                if i == 0:
                    ax.set_title(f"k = {k:g}")
        axes[0][0].set_title("image")
        for ax in axes.ravel():
            ax.set_xticks([])
            ax.set_yticks([])
        fig.tight_layout()
        figs["montage"].append((bi, _png(fig)))
    return figs


_CSS = """
:root { --fg:#1b1f24; --muted:#5b6470; --bg:#fff; --card:#f5f7fa; --line:#d9dee5;
        --accent:#2f6fb3; --hi:#fff4d6; }
@media (prefers-color-scheme: dark) {
  :root { --fg:#e6e9ee; --muted:#a3acb8; --bg:#14171b; --card:#1d2127;
          --line:#30363f; --accent:#7fb2ea; --hi:#3a3220; } }
body { font-family: system-ui, sans-serif; color: var(--fg); background: var(--bg);
       max-width: 1100px; margin: 0 auto; padding: 24px 16px; line-height: 1.5; }
h1 { font-size: 1.6em; margin-bottom: 0.2em; } h2 { margin-top: 1.8em; }
.muted { color: var(--muted); }
table { border-collapse: collapse; margin: 0.8em 0; font-size: 0.92em; }
th, td { border-bottom: 1px solid var(--line); padding: 4px 10px; text-align: right; }
th:first-child, td:first-child { text-align: left; }
tr.hi td { background: var(--hi); }
code { background: var(--card); padding: 1px 4px; border-radius: 3px; }
.params { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
          gap: 12px; }
.card { background: var(--card); border-radius: 8px; padding: 10px 14px; }
.card h3 { margin: 0 0 4px 0; font-size: 1em; color: var(--accent); }
.card p { margin: 0.3em 0; font-size: 0.92em; }
img { max-width: 100%; background: #fff; border-radius: 4px; }
.scroll { overflow-x: auto; }
"""


def _pct(v):
    if v is None or not np.isfinite(v):
        return "&gt; max" if v is not None and v == np.inf else "–"
    return f"{100 * v:.0f}%"


def render_html(info, summary, k_grid, spike_sds_um, target, census, figs, reference):
    s_list = list(summary)
    ratio = info["ratio"]
    tsd = target["spike_sd_um"]
    rows = []
    for s in s_list:
        sm = summary[s]
        ks = sm["k_suggested"]
        ats = sm.get("at_suggested", {})
        ki = None
        if np.isfinite(ks):
            ki = int(np.argmin(np.abs(np.asarray(k_grid) - ks)))
        floor = "–" if sm["k_noise_floor"] is None else f"{sm['k_noise_floor']:g}"
        cells = [
            f"<td>{s:g}</td>",
            f"<td>{sm['robust_sd']:.3g}</td>",
            f"<td>{sm['k_for_target']:.1f}</td>",
            f"<td>{floor}</td>",
            f"<td><b>{'–' if not np.isfinite(ks) else f'{ks:g}'}</b></td>",
            f"<td><code>{sm['method'] or '–'}</code></td>",
            (
                f"<td>{sm['median'] + ks * sm['robust_sd']:.4g}</td>"
                if np.isfinite(ks)
                else "<td>–</td>"
            ),
        ]
        for sd_um in spike_sds_um:
            a = ats.get(str(sd_um), {})
            cells.append(f"<td>{_pct(a.get('c50'))} / {_pct(a.get('c90'))}</td>")
        cells.append(f"<td style='text-align:left'>{sm['status']}</td>")
        cells.append(
            f"<td>{sm['objects_per_mm3'][ki]:,.0f}</td>"
            if ki is not None
            else "<td>–</td>"
        )
        hi = ' class="hi"' if s == reference.get("s") else ""
        rows.append(f"<tr{hi}>" + "".join(cells) + "</tr>")
    spike_heads = "".join(
        f"<th>sd {sd:g} µm puncta<br>50% / 90% at</th>" for sd in spike_sds_um
    )
    sugg_table = (
        "<div class='scroll'><table><tr><th>s (µm)</th><th>robust sd</th>"
        "<th>k for target</th><th>k noise floor</th><th>suggested k</th>"
        "<th>method</th><th>threshold</th>"
        + spike_heads
        + "<th>status</th><th>objects / mm³<br>(nearest k)</th></tr>"
        + "".join(rows)
        + "</table></div>"
    )

    def curve_table(key, fmt, label):
        head = (
            "<tr><th>s \\ k</th>" + "".join(f"<th>{k:g}</th>" for k in k_grid) + "</tr>"
        )
        body = ""
        for s in s_list:
            vals = summary[s][key]
            if vals is None:
                continue
            body += (
                f"<tr><td>{s:g}</td>"
                + "".join(f"<td>{fmt(v)}</td>" for v in vals)
                + "</tr>"
            )
        return f"<h3>{label}</h3><div class='scroll'><table>{head}{body}</table></div>"

    def limit_table(sd_um, which):
        head = (
            "<tr><th>s \\ k</th>" + "".join(f"<th>{k:g}</th>" for k in k_grid) + "</tr>"
        )
        body = "".join(
            f"<tr><td>{s:g}</td>"
            + "".join(f"<td>{_pct(v)}</td>" for v in summary[s]["spikes"][sd_um][which])
            + "</tr>"
            for s in s_list
        )
        n = summary[s_list[0]]["spikes"][sd_um]["n"]
        return (
            f"<h3>Contrast for {which[1:]}% detection, puncta of sd {sd_um:g} µm "
            f"<span class='muted'>({n} puncta)</span></h3>"
            f"<div class='scroll'><table>{head}{body}</table></div>"
        )

    ref_note = ""
    if reference.get("method"):
        ref_note = (
            f"<p>Your configured method <code>{reference['method']}</code> is "
            "highlighted in the tables.</p>"
        )
    census_txt = "no blobs found"
    if census.size:
        q = np.percentile(census, [25, 50, 75])
        census_txt = (
            f"median blob sd <b>{q[1]:.1f} µm</b> (IQR {q[0]:.1f}–{q[2]:.1f} µm, "
            f"{census.size} blobs)"
        )
    zooms = info["zooms_um"]
    montage_html = "".join(
        f"<h3>Block {info['blocks'][bi]['index']} "
        f"<span class='muted'>({info['blocks'][bi]['density']:,.0f} objects/mm³ at "
        f"the reference setting)</span></h3><img src='data:image/png;base64,{png}'>"
        for bi, png in figs["montage"]
    )

    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DoG tuning report</title><style>{_CSS}</style></head><body>
<h1>DoG tuning report</h1>
<p class="muted">{info['image']}<br>
stain {info['stain']} · level {info['level']} · voxel {zooms[0]:.2f} × {zooms[1]:.2f} × {zooms[2]:.2f} µm (z, y, x)
· correction <code>{info['correction_method']}</code> · {len(info['blocks'])} blocks of
{info['chunk_size']}³ voxels · k units from {info['stat_blocks']} blocks, as dogthresh</p>

<h2>What the parameters mean</h2>
<p>The DoG method is <code>dog+s&lt;s&gt;r&lt;r&gt;k&lt;k&gt;</code>; <code>p</code> is the
decimal point, so <code>dog+s2p5r4k6</code> is s = 2.5.</p>
<div class="params">
<div class="card"><h3>s — blob scale (µm)</h3><p>Sd of the inner Gaussian blur.
The DoG responds most to blobs of about this size; smaller s finds smaller, fainter
puncta and splits touching ones, larger s merges them and suppresses fine texture.</p></div>
<div class="card"><h3>r — background ratio (here {ratio:g})</h3><p>The outer Gaussian, of sd
r·s, estimates the local background that each blob is compared with. The DoG is
G(s) − G(r·s): blob minus surroundings.</p></div>
<div class="card"><h3>k — threshold (robust sd)</h3><p>A voxel is in the mask when its
DoG exceeds median + k × robust sd of the DoG in tissue (robust sd =
1.4826 × MAD). The robust sd is set mostly by real tissue texture, not camera
noise, so k is a <i>contrast</i> dial: how far above its surroundings a punctum must
be to count.</p></div>
<div class="card"><h3>Pool exclusion</h3><p>Large bright regions (stain pooling) are
found on a coarse level and removed from the mask, as in the pipeline; all counts
here exclude them.</p></div>
<div class="card"><h3>Contrast</h3><p>Peak height of a punctum relative to the local
background intensity (12 µm blur): 100% means the punctum centre is twice as bright
as its surroundings. The synthetic puncta are 3D Gaussians placed on background
tissue.</p></div>
<div class="card"><h3>Detection limit</h3><p>The contrast at which 50% or 90% of the
synthetic puncta are detected. Lower is more sensitive.</p></div>
<div class="card"><h3>Noise FDR</h3><p>Objects found in a simulated image with no
structure, only noise at the measured level, divided by objects found in the real
image. Near zero means noise is not what limits k.</p></div>
<div class="card"><h3>Scale census</h3><p>The sd of the strongest blobs in the image,
from the scale-normalised Laplacian of Gaussian. It cannot resolve blobs much smaller
than the z voxel ({zooms[0]:.1f} µm).</p></div>
</div>

<h2>Suggested settings</h2>
<p>Target: detect <b>{100 * target['recall']:.0f}%</b> of puncta of sd
<b>{tsd:g} µm</b> at <b>{100 * target['contrast']:.0f}%</b> contrast, with noise
FDR ≤ {target['max_noise_fdr']:g}. For each s, the suggested k is the largest k
(rounded down to 0.5) that meets the target. The noise floor is the lowest k on the
grid whose noise FDR is within the limit; when the target needs a k below it, the
target cannot be met at that s and no k is suggested. Larger s is less sensitive to
small puncta, so a target set for small puncta is often met only at small s.
Change the target with <code>dog_tune_target_*</code> in the config.</p>
{ref_note}
{sugg_table}
<p>Scale census: {census_txt}.</p>

<h2>Detection limits</h2>
<img src="data:image/png;base64,{figs['limits']}">
{''.join(limit_table(sd, w) for sd in spike_sds_um for w in ('c50', 'c90'))}

<h2>Object counts and noise</h2>
<img src="data:image/png;base64,{figs['counts']}">
{curve_table('objects_per_mm3', lambda v: f'{v:,.0f}', 'Objects per mm³ of tissue')}
{curve_table('volume_fraction', lambda v: f'{100 * v:.2f}%', 'Masked volume, % of tissue')}
{curve_table('noise_fdr', lambda v: f'{v:.3f}', 'Noise FDR')}

<h2>Scale census</h2>
<img src="data:image/png;base64,{figs['census']}">

<h2>Montages</h2>
<p>Maximum projection over 3 central z slices of the block core; masks in orange.
Blocks with the lowest, median and highest object density at the reference setting
(<code>{reference['label']}</code>).</p>
{montage_html}

<h2>Caveats</h2>
<ul>
<li>Synthetic puncta sit on background tissue, away from existing detections; puncta
inside dense clusters are harder to separate than these numbers suggest.</li>
<li>Contrast is measured on the corrected image. With <code>--correction_method
log</code> it is relative to log intensities and is not comparable to the other
methods.</li>
<li>k units depend on which blocks are sampled and on the image, so the same k means
a different contrast in another subject. Tune on representative subjects, then fix
s and k for the study.</li>
</ul>
</body></html>
"""


# --------------------------------------------------------------------------- #
# Snakemake entry point
# --------------------------------------------------------------------------- #


def _json_ready(obj):
    if isinstance(obj, dict):
        return {
            str(k): _json_ready(v) for k, v in obj.items() if not str(k).startswith("_")
        }
    if isinstance(obj, (list, tuple)):
        return [_json_ready(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _json_ready(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    return obj


def main():
    import dask
    import dask.array as da
    from dogthresh import load_pool_lowres, parse_dog_method, sample_dog_stats
    from zarrnii import ZarrNii

    params = snakemake.params
    ratio = float(params.ratio)
    s_grid = [float(s) for s in params.s_um]
    k_grid = [float(k) for k in params.k]
    montage_k = [float(k) for k in params.montage_k]
    spike_sds = [float(s) for s in params.spike_sd_um]
    target = dict(
        contrast=float(params.target_contrast),
        recall=float(params.target_recall),
        spike_sd_um=float(params.target_spike_sd_um),
        max_noise_fdr=float(params.max_noise_fdr),
    )
    if target["spike_sd_um"] not in spike_sds:
        spike_sds.append(target["spike_sd_um"])

    # reference setting: the first configured dog method with this ratio
    reference = dict(s=3.0, k=6.0, method=None)
    for m in params.seg_methods:
        try:
            s_ref, r_ref, k_ref = parse_dog_method(m)
        except ValueError:
            continue
        if r_ref == ratio:
            reference = dict(s=s_ref, k=k_ref, method=m)
            break
    if reference["s"] not in s_grid:
        s_grid = sorted(s_grid + [reference["s"]])
    if reference["k"] not in montage_k:
        montage_k = sorted(montage_k + [reference["k"]])
    reference["label"] = method_string(reference["s"], ratio, reference["k"])

    with dask.config.set(scheduler="threads", num_workers=snakemake.threads):
        znimg = ZarrNii.from_file(snakemake.input.corrected)
        zooms_um = [z * 1000.0 for z in znimg.get_zooms()]
        lo = load_pool_lowres(snakemake.input.corrected, params)
        chunk = int(params.chunk_size)
        img = znimg.data[0].rechunk((chunk, chunk, chunk))

        # k units exactly as dogthresh, for every s
        stats = {}
        stat_blocks = None
        for s_um in s_grid:
            stat_blocks, values = sample_dog_stats(
                img,
                lo["lowres"],
                sigma_vox(s_um, zooms_um),
                sigma_vox(s_um * ratio, zooms_um),
                int(params.stat_blocks),
                int(params.stat_samples_per_block),
            )
            _, med, sd = robust_threshold(values, 0)
            stats[s_um] = (med, sd)
            print(f"s = {s_um:g} um: DoG median {med:.4g}, robust sd {sd:.4g}")
        del values

        # analysis blocks: the first of the (randomly ordered, tissue-rich)
        # statistics blocks, with context for the largest kernel
        n_blocks = min(int(params.blocks), len(stat_blocks))
        tune_blocks = stat_blocks[:n_blocks]
        pool_struct = ellipsoid(float(params.pool_dilate_um), zooms_um)
        pool_vox = (
            [0, 0, 0] if pool_struct is None else [n // 2 for n in pool_struct.shape]
        )
        largest = max(max(s_grid) * ratio, max(params.census_scales_um))
        margin = tuple(
            int(np.ceil(_TRUNCATE * sv)) + p
            for sv, p in zip(sigma_vox(largest, zooms_um), pool_vox)
        )
        overlapped = da.overlap.overlap(
            img, depth=margin, boundary=dict.fromkeys(range(3), 0)
        )
        starts = [np.concatenate([[0], np.cumsum(c)]) for c in img.chunks]
        tasks, info_blocks = [], []
        for i, b in enumerate(tune_blocks):
            block = overlapped.blocks[b]
            start = [int(st[j]) - m for st, j, m in zip(starts, b, margin)]
            pool_up = upsample_nearest_block(lo["pool"], block.shape, start, img.shape)
            info_blocks.append(dict(index=list(map(int, b))))
            tasks.append(
                dask.delayed(analyse_block)(
                    block,
                    pool_up,
                    margin,
                    zooms_um,
                    stats,
                    ratio,
                    k_grid,
                    montage_k,
                    spike_sds,
                    int(params.spikes_per_block),
                    float(params.spike_sep_um),
                    lo["pool_thresh"],
                    float(params.pool_dilate_um),
                    [float(s) for s in params.census_scales_um],
                    int(params.census_per_block),
                    _RANDOM_STATE + i,
                )
            )
        del lo
        # each task holds a few hundred MB, so cap the parallelism
        blocks = list(
            dask.compute(
                *tasks, num_workers=min(snakemake.threads, int(params.workers))
            )
        )

    summary = summarise(blocks, stats, ratio, k_grid, spike_sds, target)
    census = np.concatenate([b["census_blob_sd"] for b in blocks])

    # montage blocks: lowest, median and highest density at the reference
    ref_s = reference["s"]
    ki = int(np.argmin(np.abs(np.asarray(k_grid) - reference["k"])))
    for ib, b in zip(info_blocks, blocks):
        r = b["s"][ref_s]
        ib["density"] = r["n_objects"][ki] / max(r["valid_mm3"], 1e-12)
        ib["noise_sd"] = b["noise_sd"]
        ib["tissue_fraction"] = b["tissue_fraction"]
    order = np.argsort([ib["density"] for ib in info_blocks])
    montage_ids = list(dict.fromkeys([order[0], order[len(order) // 2], order[-1]]))

    info = dict(
        image=str(snakemake.input.corrected),
        stain=snakemake.wildcards.stain,
        level=snakemake.wildcards.level,
        zooms_um=zooms_um,
        correction_method=params.correction_method,
        chunk_size=chunk,
        stat_blocks=int(params.stat_blocks),
        ratio=ratio,
        blocks=info_blocks,
    )
    figs = _figures(summary, k_grid, spike_sds, census, blocks, montage_ids, montage_k)
    html = render_html(
        info, summary, k_grid, spike_sds, target, census, figs, reference
    )
    with open(snakemake.output.html, "w") as f:
        f.write(html)
    with open(snakemake.output.json, "w") as f:
        json.dump(
            _json_ready(
                dict(
                    info=info,
                    target=target,
                    reference=reference,
                    k=k_grid,
                    spike_sd_um=spike_sds,
                    census_blob_sd_um=census,
                    summary={f"{s:g}": v for s, v in summary.items()},
                )
            ),
            f,
            indent=1,
        )


if __name__ == "__main__":
    main()
