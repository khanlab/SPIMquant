"""Slices and metrics of one subject in template space, for the group QC report.

Writes one JSON with the template-space slices of every stain and field fraction
map (embedded as WebP), the subject's QC metrics, and the paths of the subject's
own QC figures and reports, which the group report links to or embeds.
"""

import glob
import json
import os
import re

import nibabel as nib
import numpy as np
import pandas as pd

from groupqc_lib import (
    display_window,
    gray_u8,
    plane_axis,
    slice_indices,
    to_view,
    uri,
    value_u8,
)

# subject QC outputs linked from the report: suffix -> section, label
FIGURES = {
    "regqc.html": ("registration", "Registration QC report"),
    "histogram.png": ("intensity", "Intensity histogram"),
    "segslices.png": ("segmentation", "Segmentation slices"),
    "zprofile.png": ("segmentation", "Z-profile"),
    "objectstats.png": ("segmentation", "Object statistics"),
    "roimontage.png": ("segmentation", "ROI zoom montage"),
    "roisummary.png": ("segmentation", "ROI summary"),
    "vesselslices.png": ("vessels", "Vessel slices"),
    "vesselroimontage.png": ("vessels", "Vessel ROI montage"),
    "instancerandom.gif": ("instances", "Instances (random)"),
    "instancesorted.gif": ("instances", "Instances (by size)"),
    "colocinstancerandom.gif": ("instances", "Colocalized (random)"),
    "colocinstancesorted.gif": ("instances", "Colocalized (by size)"),
}
FIGURE_ENTITIES = ["stain", "seg", "obj", "desc"]


def _load(path):
    img = nib.as_closest_canonical(nib.load(path))
    return np.asarray(img.dataobj, dtype=np.float32).squeeze(), img


def _ncc(a, b, mask):
    a, b = a[mask].astype(np.float64), b[mask].astype(np.float64)
    if a.size < 2 or a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def _jacobian(xfm_path, mask):
    # regqc_lib needs matplotlib only for its colormaps; imported here so the
    # rest of this script does not depend on it
    from regqc_lib import JACOBIAN_STRIDE, jacobian_determinant

    det = jacobian_determinant(nib.load(xfm_path))
    sub = mask[::JACOBIAN_STRIDE, ::JACOBIAN_STRIDE, ::JACOBIAN_STRIDE]
    if det.shape != sub.shape:
        return {}
    values = det[sub]
    if values.size == 0:
        return {}
    return {
        "reg_jac_p1": float(np.percentile(values, 1)),
        "reg_jac_p99": float(np.percentile(values, 99)),
        "reg_jac_folded": float((values <= 0).mean() * 100),
    }


def _table_metrics(path, stains):
    """Whole-brain totals of each stain from a merged segstats table."""
    df = pd.read_csv(path, sep="\t")
    if "index" in df:
        df = df[df["index"] != 0]
    out = {}
    for stain in stains:
        count, volume, ff = (f"{stain}+{m}" for m in ("count", "volume", "fieldfrac"))
        if volume not in df:
            continue
        vol = pd.to_numeric(df[volume], errors="coerce")
        total = vol.sum()
        if not total > 0:
            continue
        if count in df:
            n = pd.to_numeric(df[count], errors="coerce").sum()
            out[f"{stain}:count"] = float(n)
            out[f"{stain}:density"] = float(n / total)
        if ff in df:
            weighted = (pd.to_numeric(df[ff], errors="coerce") * vol).sum()
            out[f"{stain}:fieldfrac"] = float(weighted / total)
    return out


def _figures(prefix, regqc):
    """The subject's QC figures and reports that exist, with what they show.

    prefix: path prefix of the subject's files in its qc folder; regqc: glob
    patterns of its registration QC reports.
    """
    found = []
    reports = sorted(p for pattern in regqc for p in glob.glob(pattern))
    for path in sorted(glob.glob(prefix + "_*")) + reports:
        name = os.path.basename(path)
        suffix = name.rsplit("_", 1)[-1]
        if suffix not in FIGURES or not os.path.exists(path):
            continue
        section, label = FIGURES[suffix]
        ents = dict(re.findall(r"(?:^|_)([a-zA-Z]+)-([^_]+)", name))
        desc = ents.get("desc", "")
        detail = [ents[e] for e in ("stain", "seg") if e in ents]
        if desc:
            detail.append(desc.replace("nomask", ""))
        if desc.endswith("nomask"):
            label += ", no mask"
        if suffix == "regqc.html" and desc == "tune":
            label = "Registration tuning report"
        found.append(
            {
                "section": section,
                "label": label,
                "detail": " · ".join(detail),
                "stain": ents.get("stain"),
                "path": os.path.abspath(path),
            }
        )
    return found


def build_subject(
    label,
    entities,
    template_mask,
    template_anat,
    stains,
    warped,
    seg_layers,
    fieldfrac,
    n_slices,
    xfm=None,
    segstats=None,
    figure_prefix=None,
    regqc=(),
    quality=80,
):
    mask_data, _ = _load(template_mask)
    mask = mask_data > 0
    slices = slice_indices(mask, n_slices)
    template, _ = _load(template_anat)

    layers, metrics = {}, {}
    for i, (stain, path) in enumerate(zip(stains, warped)):
        data, _ = _load(path)
        if data.shape != mask.shape:
            raise ValueError(f"{path} is not on the template grid {mask.shape}")
        window = display_window(data, mask)
        layers[f"stain:{stain}"] = {
            "kind": "gray",
            "slices": {
                plane: [
                    uri(gray_u8(to_view(data, plane_axis(plane), k), window), quality)
                    for k in idx
                ]
                for plane, idx in slices.items()
            },
        }
        inside = data[mask]
        nonzero = inside[inside > 0]
        if nonzero.size:
            metrics[f"{stain}:median"] = float(np.median(nonzero))
            metrics[f"{stain}:p99"] = float(np.percentile(nonzero, 99))
        if i == 0:
            metrics["reg_ncc_template"] = _ncc(data, template, mask)
            metrics["reg_coverage"] = float((inside > 0).mean() * 100)
        del data

    for layer, path in zip(seg_layers, fieldfrac):
        data, _ = _load(path)
        inside = data[mask]
        scale = (
            float(np.percentile(inside[inside > 0], 99.5))
            if (inside > 0).any()
            else 1.0
        )
        scale = max(scale, 1e-3)
        layers[f"fieldfrac:{layer['stain']}:{layer['desc']}"] = {
            "kind": "value",
            "stain": layer["stain"],
            "desc": layer["desc"],
            "scale": scale,
            "slices": {
                plane: [
                    uri(value_u8(to_view(data, plane_axis(plane), k), scale), 90)
                    for k in idx
                ]
                for plane, idx in slices.items()
            },
        }
        metrics[f"{layer['stain']}:{layer['desc']}:fieldfrac_mean"] = float(
            inside.mean()
        )
        del data

    if xfm:
        metrics.update(_jacobian(xfm, mask))
    if segstats:
        metrics.update(_table_metrics(segstats, [l["stain"] for l in seg_layers]))

    return {
        "label": label,
        "entities": entities,
        "slices": slices,
        "layers": layers,
        "metrics": {k: v for k, v in metrics.items() if v is not None},
        "figures": _figures(figure_prefix, regqc) if figure_prefix else [],
    }


if "snakemake" in globals():
    p, inp = snakemake.params, snakemake.input
    result = build_subject(
        label=p.label,
        entities=p.entities,
        template_mask=inp.template_mask,
        template_anat=inp.template_anat,
        stains=p.stains,
        warped=list(inp.warped),
        seg_layers=p.seg_layers,
        fieldfrac=list(inp.fieldfrac),
        n_slices=p.n_slices,
        xfm=inp.get("xfm") or None,
        segstats=inp.get("segstats") or None,
        figure_prefix=p.figure_prefix,
        regqc=p.regqc,
        quality=p.quality,
    )
    with open(snakemake.output.json, "w") as f:
        json.dump(result, f)
