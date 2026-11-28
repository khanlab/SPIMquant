"""Group QC report: one self-contained HTML page on every subject of a run.

Takes the per-subject JSONs of groupqc_subject.py and adds what needs the whole
group: the template slices and atlas outlines, each subject's correlation with
the mean of the others, a color scale shared by all subjects for each field
fraction map, and outlier flags (robust z-scores against the group).

The page has four views: one subject at a time (with its registration,
stains, segmentation and its own QC figures, stepping through subjects with the
arrow keys), a montage of all subjects at one slice, an animated flip through
the subjects, and a sortable table of the metrics.
"""

import base64
import datetime
import html
import json
import os
import re
from io import BytesIO

import nibabel as nib
import numpy as np
import pandas as pd
from PIL import Image

import groupqc_lib
from groupqc_lib import (
    METRICS,
    PLANES,
    SEG_METRICS,
    STAIN_METRICS,
    decode_uri,
    display_window,
    gray_u8,
    outline_uri,
    plane_axis,
    slice_indices,
    to_view,
    uri,
)

MIN_SUBJECTS_FOR_OUTLIERS = 5
MIN_RELATIVE_SPREAD = 0.01
EMBED_MAX_PX = 1600
EMBED_MAX_GIF_MB = 4
ASSETS = os.path.dirname(os.path.abspath(groupqc_lib.__file__))


def _load(path):
    img = nib.as_closest_canonical(nib.load(path))
    return np.asarray(img.dataobj).squeeze(), img


def _metric_info(key):
    """Label, section, unit, direction, help and stain of a metric key."""
    if key in METRICS:
        label, section, unit, direction, text = METRICS[key]
        return dict(label=label, section=section, unit=unit, dir=direction, help=text)
    parts = key.split(":")
    if len(parts) == 2 and parts[1] in STAIN_METRICS:
        label, unit, direction, text = STAIN_METRICS[parts[1]]
        return dict(
            label=f"{parts[0]} {label.lower()}",
            section="intensity",
            unit=unit,
            dir=direction,
            help=text,
            stain=parts[0],
        )
    if len(parts) == 2 and parts[1] in SEG_METRICS:
        label, unit, direction, text = SEG_METRICS[parts[1]]
        return dict(
            label=f"{parts[0]} {label.lower()}",
            section="segmentation",
            unit=unit,
            dir=direction,
            help=text,
            stain=parts[0],
        )
    if len(parts) == 3 and parts[2] in SEG_METRICS:
        label, unit, direction, text = SEG_METRICS[parts[2]]
        return dict(
            label=f"{parts[0]} {label.lower()}",
            section="segmentation",
            unit=unit,
            dir=direction,
            help=text,
            stain=parts[0],
            desc=parts[1],
        )
    return None


def _robust(values):
    values = np.asarray([v for v in values if v is not None], dtype=float)
    if values.size == 0:
        return None
    median = float(np.median(values))
    spread = float(1.4826 * np.median(np.abs(values - median)))
    # when most subjects share a value (e.g. 100% coverage), the MAD is ~0 and
    # every small difference would be an outlier: spread is at least 1% of the median
    spread = max(spread, MIN_RELATIVE_SPREAD * abs(median))
    return {"median": median, "spread": spread, "n": int(values.size)}


def _ncc_to_group(subjects, layer_key, masks):
    """Correlation of each subject's slices with the mean of the others."""
    stacks = {}
    for i, s in enumerate(subjects):
        layer = s["layers"].get(layer_key)
        if layer is None:
            continue
        stacks[i] = np.concatenate(
            [
                decode_uri(u).astype(np.float32)[masks[plane][k]]
                for plane, uris in layer["slices"].items()
                for k, u in enumerate(uris)
            ]
        )
    if len(stacks) < 3:
        return {}
    total = np.sum(list(stacks.values()), axis=0)
    out = {}
    for i, v in stacks.items():
        others = (total - v) / (len(stacks) - 1)
        if v.std() > 0 and others.std() > 0:
            out[i] = float(np.corrcoef(v, others)[0, 1])
    return out


def _figure_src(path, html_dir, embed):
    if not embed or path.endswith(".html"):
        return os.path.relpath(path, html_dir)
    if path.endswith(".gif"):
        if os.path.getsize(path) > EMBED_MAX_GIF_MB * 1e6:
            return os.path.relpath(path, html_dir)
        with open(path, "rb") as f:
            return "data:image/gif;base64," + base64.b64encode(f.read()).decode()
    image = Image.open(path)
    image.thumbnail((EMBED_MAX_PX, EMBED_MAX_PX))
    buffer = BytesIO()
    image.convert("RGB").save(buffer, format="WEBP", quality=75, method=4)
    return "data:image/webp;base64," + base64.b64encode(buffer.getvalue()).decode()


def _anchor(label):
    return re.sub(r"[^A-Za-z0-9-]+", "-", label).strip("-")


def _covariates(participants_tsv, subjects):
    if not participants_tsv or not os.path.exists(participants_tsv):
        return [{} for _ in subjects]
    df = pd.read_csv(participants_tsv, sep="\t", dtype=str).fillna("")
    df = df.set_index("participant_id") if "participant_id" in df else df
    out = []
    for s in subjects:
        pid = "sub-" + str(s["entities"].get("subject", ""))
        if pid in df.index:
            row = df.loc[pid]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            out.append({k: v for k, v in row.items() if v not in ("", "n/a")})
        else:
            out.append({})
    return out


def build_data(
    subject_jsons,
    template_anat,
    template_mask,
    outline_dseg,
    html_path,
    title,
    template_name,
    participants_tsv=None,
    group_by=(),
    outlier_z=3.0,
    embed_figures=False,
    reg_method="",
    version="",
    quality=80,
):
    subjects = []
    for path in subject_jsons:
        with open(path) as f:
            subjects.append(json.load(f))
    if not subjects:
        raise ValueError("No subjects to report on")
    subjects.sort(key=lambda s: s["label"])
    html_dir = os.path.dirname(os.path.abspath(html_path))

    # slices are the same for every subject; the first sets the plan
    plan = subjects[0]["slices"]
    for s in subjects:
        if s["slices"] != plan:
            raise ValueError(
                f"{s['label']} has other slice positions than {subjects[0]['label']}: "
                "were they made with the same template and settings?"
            )

    mask, mask_img = _load(template_mask)
    mask = mask > 0
    anat, _ = _load(template_anat)
    anat = anat.astype(np.float32)
    window = display_window(anat, mask)
    zooms = mask_img.header.get_zooms()[:3]
    labels = _load(outline_dseg)[0].astype(np.int32) if outline_dseg else None
    if labels is not None and labels.shape != mask.shape:
        from nibabel.processing import resample_from_to

        labels = np.asarray(
            resample_from_to(
                nib.as_closest_canonical(nib.load(outline_dseg)), mask_img, order=0
            ).dataobj
        ).astype(np.int32)

    planes, tpl_anat, tpl_outlines, view_masks = {}, {}, {}, {}
    for plane, indices in plan.items():
        axis = plane_axis(plane)
        direction = next(d for p, a, d in PLANES if p == plane)
        a, b = [i for i in range(3) if i != axis]
        width = mask.shape[a] * zooms[a]
        height = mask.shape[b] * zooms[b]
        origin, spacing = mask_img.affine[axis, 3], mask_img.affine[axis, axis]
        planes[plane] = {
            "aspect": f"{width:.4f} / {height:.4f}",
            "positions": [
                f"{direction} {origin + k * spacing:+.2f} mm" for k in indices
            ],
        }
        tpl_anat[plane] = [
            uri(gray_u8(to_view(anat, axis, k), window), quality) for k in indices
        ]
        view_masks[plane] = [to_view(mask, axis, k) for k in indices]
        if labels is not None:
            tpl_outlines[plane] = [
                outline_uri(to_view(labels, axis, k)) for k in indices
            ]

    # layers in the order of the first subject that has them: stains, then maps
    layers = {}
    for s in subjects:
        for key, layer in s["layers"].items():
            if key in layers:
                continue
            kind, stain = layer["kind"], key.split(":")[1]
            if kind == "gray":
                layers[key] = dict(
                    key=key, kind=kind, stain=stain, label=stain, group="Stains"
                )
            else:
                layers[key] = dict(
                    key=key,
                    kind=kind,
                    stain=stain,
                    desc=layer.get("desc", ""),
                    label=f"{stain} field fraction",
                    group="Field fraction",
                    unit="%",
                )
    gray_first = sorted(layers.values(), key=lambda l: l["kind"] != "gray")
    reg_key = next((l["key"] for l in gray_first if l["kind"] == "gray"), None)

    if reg_key:
        for i, ncc in _ncc_to_group(subjects, reg_key, view_masks).items():
            subjects[i]["metrics"]["reg_ncc_group"] = ncc

    shared_scale = {
        key: float(
            np.median(
                [s["layers"][key]["scale"] for s in subjects if key in s["layers"]]
            )
        )
        for key, l in layers.items()
        if l["kind"] == "value"
    }

    # metrics: registration first, then by stain
    keys = list(dict.fromkeys(k for s in subjects for k in s["metrics"]))
    order = {"registration": 0, "intensity": 1, "segmentation": 2}
    metrics = []
    for key in keys:
        info = _metric_info(key)
        if info is None:
            continue
        metrics.append(dict(key=key, **info))
    reg_order = list(METRICS)
    stain_order = [l["stain"] for l in gray_first]
    metrics.sort(
        key=lambda m: (
            order.get(m["section"], 9),
            reg_order.index(m["key"]) if m["key"] in METRICS else 0,
            stain_order.index(m["stain"]) if m.get("stain") in stain_order else 99,
        )
    )

    stats, out_subjects = {}, []
    for m in metrics:
        st = _robust(s["metrics"].get(m["key"]) for s in subjects)
        if st:
            stats[m["key"]] = st
    covs = _covariates(participants_tsv, subjects)
    for s, cov in zip(subjects, covs):
        z, flags = {}, {}
        for m in metrics:
            v, st = s["metrics"].get(m["key"]), stats.get(m["key"])
            if v is None or not st or st["spread"] == 0:
                continue
            zv = (v - st["median"]) / st["spread"]
            z[m["key"]] = round(zv, 2)
            if st["n"] >= MIN_SUBJECTS_FOR_OUTLIERS and (
                (m["dir"] in ("high", "both") and zv > outlier_z)
                or (m["dir"] in ("low", "both") and zv < -outlier_z)
            ):
                flags[m["key"]] = f"{'high' if zv > 0 else 'low'} (z = {zv:+.1f})"
        figures = [
            dict(
                section=f["section"],
                label=f["label"],
                detail=f["detail"],
                path=os.path.basename(f["path"]),
                src=_figure_src(f["path"], html_dir, embed_figures),
            )
            for f in s["figures"]
            if os.path.exists(f["path"])
        ]
        out_subjects.append(
            dict(
                label=s["label"],
                anchor=_anchor(s["label"]),
                cov=cov,
                metrics=s["metrics"],
                z=z,
                flags=flags,
                layers={
                    k: {"scale": l.get("scale"), "slices": l["slices"]}
                    for k, l in s["layers"].items()
                },
                figures=figures,
            )
        )

    present = [c for c in group_by if any(c in cov for cov in covs)]
    stains = [l["stain"] for l in gray_first if l["kind"] == "gray"]
    segs = sorted(
        {f"{l['stain']} {l['desc']}" for l in layers.values() if l["kind"] == "value"}
    )
    return {
        "title": title,
        "template_name": template_name,
        "facts": [
            ("Subjects", str(len(subjects))),
            ("Template", template_name),
            ("Stains", ", ".join(stains)),
            *([("Registration", reg_method)] if reg_method else []),
            *([("Segmentation", "; ".join(segs))] if segs else []),
            *([("SPIMquant", version)] if version else []),
            ("Written", datetime.date.today().isoformat()),
        ],
        "planes": planes,
        "template": {"anat": tpl_anat, "outlines": tpl_outlines or None},
        "layers": gray_first,
        "metrics": metrics,
        "metric_stats": stats,
        "shared_scale": shared_scale,
        "outlier_z": outlier_z,
        "group_by": present,
        "embedded": bool(embed_figures),
        "subjects": out_subjects,
    }


def metrics_tsv(data, path):
    rows = []
    for s in data["subjects"]:
        row = {
            "subject": s["label"],
            **{c: s["cov"].get(c, "") for c in data["group_by"]},
        }
        row["n_flags"] = len(s["flags"])
        row["flags"] = "; ".join(f"{k} {why}" for k, why in s["flags"].items())
        for m in data["metrics"]:
            row[m["key"]] = s["metrics"].get(m["key"])
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, na_rep="n/a")


def render(data, html_path):
    with open(os.path.join(ASSETS, "groupqc_report.css")) as f:
        css = f.read()
    with open(os.path.join(ASSETS, "groupqc_report.js")) as f:
        js = f.read()
    facts = "".join(
        f"<span>{html.escape(k)} <b>{html.escape(v)}</b></span>"
        for k, v in data["facts"]
    )
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    title = html.escape(data["title"])
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="claude:viewport-fit-sides" content="cover">
<title>{title} group QC</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&family=Space+Grotesk:wght@500;600&display=swap">
<style>{css}</style>
</head>
<body>
<header class="topbar">
  <div class="wrap">
    <div class="brand">
      <p class="eyebrow">SPIMquant group QC</p>
      <h1>{title}</h1>
    </div>
    <nav class="seg" id="tabs" aria-label="Views">
      <button type="button" data-view="subject">Subjects</button>
      <button type="button" data-view="montage">Montage</button>
      <button type="button" data-view="flip">Flip</button>
      <button type="button" data-view="metrics">Metrics</button>
    </nav>
    <button class="btn" type="button" id="theme" aria-label="Switch light or dark theme" title="Light or dark theme"><svg viewBox="0 0 16 16" aria-hidden="true"><circle class="ring" cx="8" cy="8" r="6.5"/><path class="fill" d="M8 1.5a6.5 6.5 0 0 0 0 13z"/></svg></button>
  </div>
  <div class="wrap facts">{facts}</div>
</header>
<main class="wrap">
  <div class="view" id="view-subject">
    <div class="subject-layout">
      <aside class="sidebar" aria-label="Subjects">
        <div class="search"><input type="search" id="side-search" placeholder="Filter subjects" aria-label="Filter subjects"></div>
        <ol id="side-list"></ol>
      </aside>
      <div id="subject-main" style="min-width:0"></div>
    </div>
  </div>
  <div class="view" id="view-montage" hidden>
    <div class="panel-controls" id="montage-controls"></div>
    <div class="montage" id="montage"></div>
    <p class="hint">Every subject at the same template slice. Tiles with a red border have flagged metrics; hover over a tile to see which.</p>
  </div>
  <div class="view" id="view-flip" hidden>
    <div class="panel-controls" id="flip-controls"></div>
    <div class="flip-stage" id="flip-stage"></div>
    <p class="hint">Flipping through aligned subjects makes misregistration, missing tissue and odd intensity stand out as jumps. Red ticks under the slider mark subjects with flagged metrics.</p>
  </div>
  <div class="view" id="view-metrics" hidden>
    <div class="dists" id="dists"></div>
    <div class="table-wrap"><table id="metrics-table"></table></div>
    <p class="hint">Flags compare each subject with the group: a robust z-score (median and MAD) beyond ±{data["outlier_z"]:g} in the direction that signals a problem. They point at what to look at; they do not decide whether a subject is usable. Click a column to sort and a row to open the subject.</p>
  </div>
</main>
<script type="application/json" id="qc-data">{payload}</script>
<script>{js}</script>
</body>
</html>
"""
    with open(html_path, "w") as f:
        f.write(page)


if "snakemake" in globals():
    p = snakemake.params
    inp = snakemake.input
    data = build_data(
        subject_jsons=list(inp.subjects),
        template_anat=inp.template_anat,
        template_mask=inp.template_mask,
        outline_dseg=inp.get("outline_dseg") or None,
        html_path=snakemake.output.html,
        title=p.title,
        template_name=p.template,
        participants_tsv=inp.get("participants") or None,
        group_by=p.group_by,
        outlier_z=p.outlier_z,
        embed_figures=p.embed_figures,
        reg_method=p.reg_method,
        version=p.version,
        quality=p.quality,
    )
    metrics_tsv(data, snakemake.output.tsv)
    render(data, snakemake.output.html)
