"""Landmark target registration error (TRE) for SPIMquant transforms.

Points are warped with antsApplyTransformsToPoints (ITK/LPS convention). Points
move opposite to images: to take points from space A to space B, apply the
transform that resamples B-space images into A (the B->A, "from-B_to-A" xfm),
or the fixed-point inverse of the A->B one.

Landmark files are 3D Slicer markups (.mrk.json / *_landmarks.json), LPS, mm.
A file may hold one or many control points; labels are matched across files by
name with a leading "Template_" and a trailing "-<n>" stripped (files are never
modified).

Run from Snakemake (snakemake object present) or from the command line:

  python landmark_tre.py --src <landmarks> --dst <landmarks> \
      --xfm <composite.nii.gz>[,1] [--xfm ...] | --fixedpoint-of <fwd.nii.gz> \
      --out out.tsv [--meta subject=X path=spim_to_template method=greedy stain=YoPro]
"""

import argparse
import csv
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

_LABEL_RE = re.compile(r"^(?:Template_)?(?:SPIM_|MRI_)?(?P<name>.+?)(?:-\d+)?$")


def landmark_name(label):
    """Template_CC_genu-1 -> CC_genu; Template_CC_genu -> CC_genu;
    Template_SPIM_striatum_L_dorsomedial -> striatum_L_dorsomedial (modality prefix)."""
    return _LABEL_RE.match(label.strip()).group("name")


def load_landmarks(path):
    """Read landmarks from a markups file or a folder of them.

    Returns {name: np.array([x, y, z])} in LPS mm. Points flagged as undefined
    (positionStatus != defined) are skipped. Duplicate names raise.
    """
    files = (
        sorted(glob.glob(os.path.join(path, "*.json")))
        if os.path.isdir(path)
        else [path]
    )
    out = {}
    for f in files:
        with open(f) as fh:
            js = json.load(fh)
        for mk in js.get("markups", []):
            cs = mk.get("coordinateSystem", "LPS")
            for cp in mk.get("controlPoints", []):
                if cp.get("positionStatus", "defined") != "defined":
                    continue
                xyz = np.asarray(cp["position"], float)
                if cs == "RAS":
                    xyz = xyz * np.array([-1.0, -1.0, 1.0])
                elif cs != "LPS":
                    raise ValueError(f"{f}: unsupported coordinateSystem {cs}")
                name = landmark_name(cp.get("label", ""))
                if name in out:
                    raise ValueError(f"duplicate landmark {name} in {path}")
                out[name] = xyz
    return out


def _ants_points_bin():
    exe = shutil.which("antsApplyTransformsToPoints")
    if exe is None:
        raise RuntimeError("antsApplyTransformsToPoints not found on PATH")
    return exe


_SITK_CACHE = {}


def _parse_xfm(t):
    """'[f,1]' / 'f,1' -> (f, True); 'f' -> (f, False)."""
    t = t.strip()
    if t.startswith("[") and t.endswith("]"):
        t = t[1:-1]
    if t.endswith(",1"):
        return t[:-2], True
    if t.endswith(",0"):
        return t[:-2], False
    return t, False


def _sitk_transform(path, invert):
    import SimpleITK as sitk

    key = (path, invert)
    if key not in _SITK_CACHE:
        if path.endswith((".nii", ".nii.gz")):
            if invert:
                raise ValueError(f"cannot invert a displacement field directly: {path}")
            img = sitk.ReadImage(path, sitk.sitkVectorFloat64)
            tx = sitk.DisplacementFieldTransform(img)
        else:
            tx = sitk.ReadTransform(path)
            if invert:
                tx = tx.GetInverse()
        _SITK_CACHE[key] = tx
    return _SITK_CACHE[key]


def field_coverage(path, points):
    """True for each LPS point inside a displacement field's grid (header only).

    Outside the grid ITK uses zero displacement, which would silently give a
    wrong position, so callers flag those points."""
    import SimpleITK as sitk

    r = sitk.ImageFileReader()
    r.SetFileName(path)
    r.ReadImageInformation()
    size = np.array(r.GetSize()[:3])
    origin, spacing = np.array(r.GetOrigin()[:3]), np.array(r.GetSpacing()[:3])
    direction = np.array(r.GetDirection()).reshape(len(r.GetSize()), -1)[:3, :3]
    idx = (np.linalg.inv(direction @ np.diag(spacing)) @ (np.atleast_2d(points) - origin).T).T
    return np.all((idx >= -0.5) & (idx <= size - 0.5), axis=1)


def apply_transforms(points, transforms, backend="sitk"):
    """Warp an (N,3) LPS array with ANTs-style transform list.

    Order follows antsApplyTransforms -t: the LAST transform is applied to the
    points first. '[f,1]' or 'f,1' inverts an affine. backend='sitk' (cached,
    fast) or 'ants' (antsApplyTransformsToPoints CLI, reference).
    """
    if backend == "sitk":
        pts = np.atleast_2d(np.asarray(points, float))
        txs = [_sitk_transform(*_parse_xfm(t)) for t in transforms]
        out = np.empty_like(pts)
        for i, p in enumerate(pts):
            q = tuple(p)
            for tx in reversed(txs):
                q = tx.TransformPoint(q)
            out[i] = q
        return out
    return _apply_transforms_ants(points, transforms)


def _apply_transforms_ants(points, transforms):
    pts = np.atleast_2d(np.asarray(points, float))
    with tempfile.TemporaryDirectory() as td:
        fin, fout = os.path.join(td, "in.csv"), os.path.join(td, "out.csv")
        np.savetxt(
            fin,
            np.c_[pts, np.zeros(len(pts))],
            delimiter=",",
            header="x,y,z,t",
            comments="",
        )
        cmd = [_ants_points_bin(), "-d", "3", "-i", fin, "-o", fout]
        for t in transforms:
            f, inv = _parse_xfm(t)
            cmd += ["-t", f"[{f},1]" if inv else f]
        subprocess.run(cmd, check=True, capture_output=True)
        return np.loadtxt(fout, delimiter=",", skiprows=1, ndmin=2)[:, :3]


def fixedpoint_inverse(points, fwd_transforms, init_transforms=None, iters=50, tol_mm=1e-4):
    """Points through the inverse of fwd_transforms by fixed-point iteration.

    fwd maps target-space points y to source-space points x (it is the transform
    that resamples source images into target space). We solve fwd(y) = x for y,
    starting from init_transforms(x) if given, else x:  y <- y - (fwd(y) - x).
    Returns (y, max_residual_mm, n_iter).
    """
    x = np.atleast_2d(np.asarray(points, float))
    y = apply_transforms(x, init_transforms) if init_transforms else x.copy()
    res = np.inf
    for it in range(1, iters + 1):
        r = apply_transforms(y, fwd_transforms) - x
        y = y - r
        res = float(np.abs(r).max())
        if res < tol_mm:
            break
    return y, res, it


def tre_rows(src, dst, warped, meta):
    """One row per landmark present in both src and dst."""
    rows = []
    for name, w in warped.items():
        if name not in dst:
            continue
        d = dst[name]
        rows.append(
            dict(
                **meta,
                landmark=name,
                tre_um=round(float(np.linalg.norm(w - d)) * 1000.0, 2),
                src_x=src[name][0], src_y=src[name][1], src_z=src[name][2],
                warped_x=w[0], warped_y=w[1], warped_z=w[2],
                dst_x=d[0], dst_y=d[1], dst_z=d[2],
            )
        )
    return rows


def compute(src_file, dst_file, transforms=None, fixedpoint_of=None, init=None, meta=None):
    src, dst = load_landmarks(src_file), load_landmarks(dst_file)
    names = [n for n in src if n in dst]
    if not names:
        raise ValueError(f"no landmark names shared by {src_file} and {dst_file}")
    X = np.array([src[n] for n in names])
    meta = dict(meta or {})
    if fixedpoint_of:
        Y, res, nit = fixedpoint_inverse(X, fixedpoint_of, init_transforms=init)
        meta.update(inversion="fixedpoint", inv_residual_um=round(res * 1000, 3), inv_iters=nit)
    else:
        Y = apply_transforms(X, transforms)
        meta.update(inversion="none", inv_residual_um="", inv_iters="")
    warped = dict(zip(names, Y))
    rows = tre_rows(src, dst, warped, meta)
    missing = sorted(set(src) ^ set(dst))
    return rows, missing


def write_tsv(rows, path, fieldnames=None):
    """Rows as a TSV. With no rows, writes the header of fieldnames (or an empty
    file without them), so readers see an empty table rather than no columns."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as f:
        fields = list(rows[0]) if rows else list(fieldnames or [])
        if not fields:
            return
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _cli():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True, help="landmarks in the moving/source space")
    ap.add_argument("--dst", required=True, help="landmarks in the target space")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--xfm", action="append", help="ANTs transform for points (repeat; ',1' = invert)")
    g.add_argument("--fixedpoint-of", action="append", help="forward transform(s) to invert by fixed-point")
    ap.add_argument("--init", action="append", help="initial transform(s) for fixed-point inversion")
    ap.add_argument("--out", required=True)
    ap.add_argument("--meta", nargs="*", default=[], help="key=value columns to add")
    a = ap.parse_args()
    meta = dict(kv.split("=", 1) for kv in a.meta)
    rows, missing = compute(a.src, a.dst, a.xfm, a.fixedpoint_of, a.init, meta)
    write_tsv(rows, a.out)
    tre = np.array([r["tre_um"] for r in rows])
    print(f"{a.out}: n={len(tre)} median={np.median(tre):.1f} um "
          f"IQR={np.percentile(tre,25):.1f}-{np.percentile(tre,75):.1f}"
          + (f"  (unmatched: {', '.join(missing)})" if missing else ""))


if __name__ == "__main__":
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        rows, _ = compute(
            sm.input.src, sm.input.dst,
            transforms=sm.params.get("xfm"),
            fixedpoint_of=sm.params.get("fixedpoint_of"),
            init=sm.params.get("init"),
            meta=sm.params.get("meta"),
        )
        write_tsv(rows, sm.output.tsv)
    else:
        _cli()
