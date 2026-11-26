"""TRE of one scan's manual landmarks for every stain x method x path (landmark_tre rule).

Inputs (named): spim_landmarks, template_landmarks, [mri_landmarks],
fwd__<stain>__<method> (from-subject_to-<tpl>, template grid),
inv__<stain>__<method> (from-<tpl>_to-subject, subject grid),
[mri2spim__<stain>] (from-SPIM_to-<mri>, MRI grid: MRI points -> SPIM),
[spim2mri__<stain>] (from-<mri>_to-SPIM, SPIM grid: SPIM points -> MRI).
"""

import sys

import numpy as np

sm = snakemake  # noqa: F821
sys.path.insert(0, sm.params.scriptdir)
import landmark_tre as LT  # noqa: E402

inp = sm.input
P = sm.params
spim = LT.load_landmarks(inp.spim_landmarks)
tpl = LT.load_landmarks(inp.template_landmarks)
mri = LT.load_landmarks(inp.mri_landmarks) if "mri_landmarks" in inp.keys() else None
scan = {k: v for k, v in sm.wildcards.items()}

rows = []
COLUMNS = list(scan) + [
    "path",
    "method",
    "stain",
    "inversion",
    "landmark",
    "tre_um",
    "inside_field",
    "inv_residual_um",
    "src_x",
    "src_y",
    "src_z",
    "warped_x",
    "warped_y",
    "warped_z",
    "dst_x",
    "dst_y",
    "dst_z",
]


def add(path, method, stain, inversion, src, dst, warped, inside, resid=None):
    if not len(warped["names"]):
        return
    for i, name in enumerate(warped["names"]):
        if name not in dst:
            continue
        w, d = warped["pts"][i], dst[name]
        rows.append(
            dict(
                **scan,
                path=path,
                method=method,
                stain=stain,
                inversion=inversion,
                landmark=name,
                tre_um=round(float(np.linalg.norm(w - d)) * 1000, 2),
                inside_field=bool(inside[i]),
                inv_residual_um=(
                    "" if resid is None else round(float(resid[i]) * 1000, 4)
                ),
                src_x=src[name][0],
                src_y=src[name][1],
                src_z=src[name][2],
                warped_x=w[0],
                warped_y=w[1],
                warped_z=w[2],
                dst_x=d[0],
                dst_y=d[1],
                dst_z=d[2],
            )
        )


def pts_of(d, names):
    return np.array([d[n] for n in names]).reshape(-1, 3)


def coverage(f, X):
    return LT.field_coverage(f, X) if len(X) else np.zeros(0, bool)


def warp(X, fs):
    return LT.apply_transforms(X, fs) if len(X) else X.copy()


def fixedpoint(X, fwd, inv):
    """Fixed-point inverse of fwd from the shipped inverse; per-point residual (mm)."""
    if not len(X):
        return X.copy(), np.zeros(0)
    Y = LT.apply_transforms(X, [inv])
    for _ in range(50):
        R = LT.apply_transforms(Y, [fwd]) - X
        Y = Y - R
        if np.abs(R).max() < 1e-5:
            break
    resid = np.linalg.norm(LT.apply_transforms(Y, [fwd]) - X, axis=1)
    return Y, resid


s_names = [n for n in spim if n in tpl]
t_names = [n for n in tpl if n in spim]
m_names = [n for n in mri if n in spim] if mri else []
if not s_names:
    print(
        f"WARNING {scan}: no landmark name is shared by the SPIM and template landmarks; "
        "see landmarks/group/landmark_name_check.tsv",
        file=sys.stderr,
    )

for stain in P.stains:
    mri_in_spim = None
    if mri:
        # MRI points -> SPIM (from-SPIM_to-MRI, defined on the MRI grid)
        f = inp[f"mri2spim__{stain}"]
        X = pts_of(mri, m_names)
        Y = warp(X, [f])
        add(
            "mri_to_spim",
            "greedy",
            stain,
            "none",
            mri,
            spim,
            dict(names=m_names, pts=Y),
            coverage(f, X),
        )
        mri_in_spim = (Y, coverage(f, X))
        # SPIM points -> MRI (from-MRI_to-SPIM, defined on the SPIM grid)
        f = inp[f"spim2mri__{stain}"]
        sm_names = [n for n in spim if n in mri]
        X = pts_of(spim, sm_names)
        add(
            "spim_to_mri",
            "greedy",
            stain,
            "none",
            spim,
            mri,
            dict(names=sm_names, pts=warp(X, [f])),
            coverage(f, X),
        )

    for method in P.methods:
        fwd, inv = inp[f"fwd__{stain}__{method}"], inp[f"inv__{stain}__{method}"]
        # SPIM -> template: shipped inverse, and fixed-point inverse of fwd
        X = pts_of(spim, s_names)
        inside = coverage(inv, X)
        Yship = warp(X, [inv])
        add(
            "spim_to_template",
            method,
            stain,
            "shipped",
            spim,
            tpl,
            dict(names=s_names, pts=Yship),
            inside,
        )
        Yfp, res = fixedpoint(X, fwd, inv)
        add(
            "spim_to_template",
            method,
            stain,
            "fixedpoint",
            spim,
            tpl,
            dict(names=s_names, pts=Yfp),
            inside & coverage(fwd, Yfp),
            res,
        )
        # template -> SPIM: forward composite (template grid) on template points
        T = pts_of(tpl, t_names)
        add(
            "template_to_spim",
            method,
            stain,
            "none",
            tpl,
            spim,
            dict(names=t_names, pts=warp(T, [fwd])),
            coverage(fwd, T),
        )
        # MRI -> template via SPIM (same stain for both legs)
        if mri_in_spim is not None:
            Ym, inside_m = mri_in_spim
            mt = [i for i, n in enumerate(m_names) if n in tpl]
            names = [m_names[i] for i in mt]
            Ys = Ym[mt]
            inside2 = inside_m[mt] & coverage(inv, Ys)
            add(
                "mri_to_template",
                method,
                stain,
                "shipped",
                mri,
                tpl,
                dict(names=names, pts=warp(Ys, [inv])),
                inside2,
            )
            Yfp2, res2 = fixedpoint(Ys, fwd, inv)
            add(
                "mri_to_template",
                method,
                stain,
                "fixedpoint",
                mri,
                tpl,
                dict(names=names, pts=Yfp2),
                inside2 & coverage(fwd, Yfp2),
                res2,
            )
        LT._SITK_CACHE.pop((fwd, False), None)  # template-grid fields are large

LT.write_tsv(rows, sm.output.tsv, fieldnames=COLUMNS)
