"""Build registration QC reports: self-contained HTML slice viewers.

Every report compares one or more candidates (registration stages or methods)
against a fixed image. Slices are pre-rendered and embedded, so the report needs
no server or network. A candidate is a dict with:

- name: label shown in the report
- warped: path of the moving image resampled onto the fixed image
- labels: optional path of an atlas segmentation in fixed image space
- warp: optional path of an ITK displacement field, for the Jacobian layer
- note: optional text describing the candidate, e.g. its parameters

A report can also have several stains, to choose between in the report. Each is a
dict with:

- name: label shown in the report
- fixed: path of the fixed image of the stain
- warped: for each candidate, the path of the moving image for this stain resampled
  onto the fixed image, which can be the same image for every stain, or None if that
  candidate has none, in which case it is hidden for the stain
"""

import base64
import html
import json
from io import BytesIO

import nibabel as nib
import numpy as np
from matplotlib import colormaps
from nibabel.processing import resample_from_to
from PIL import Image

N_SLICES = 7
JACOBIAN_STRIDE = 2
LOG_JACOBIAN_RANGE = 1.0
# atlas outlines are drawn at about this width, so they stay thin on coarse images
OUTLINE_WIDTH_PX = 1000

# plane name, array axis through the plane (RAS+ arrays), name of that direction
PLANES = [("Coronal", 1, "A"), ("Axial", 2, "S"), ("Sagittal", 0, "R")]


MEDIAN_GRAY = 0.4  # display gray level of an image's median voxel
GAMMA_RANGE = (0.2, 1.0)


def jacobian_determinant(warp_img, stride=JACOBIAN_STRIDE):
    """Jacobian determinant of an ITK displacement field (LPS mm displacements).

    Returns the determinant, computed every `stride` voxels to limit memory use.
    """
    disp = np.squeeze(np.asarray(warp_img.dataobj, dtype=np.float32))
    if disp.ndim != 4 or disp.shape[-1] != 3:
        raise ValueError(f"Expected a 3D vector field, got shape {disp.shape}")
    disp = disp[::stride, ::stride, ::stride]
    # voxel to LPS coordinates, in which ITK stores displacements
    voxel_to_lps = np.diag([-1.0, -1.0, 1.0]) @ warp_img.affine[:3, :3] * stride
    lps_to_voxel = np.linalg.inv(voxel_to_lps)
    jac = np.empty(disp.shape[:3] + (3, 3), dtype=np.float32)
    for component in range(3):
        grads = np.gradient(disp[..., component])
        for axis in range(3):
            jac[..., component, axis] = grads[axis]
    jac = jac @ lps_to_voxel.astype(np.float32)
    jac[..., [0, 1, 2], [0, 1, 2]] += 1
    det = np.linalg.det(jac)
    # same orientation as the images, which are shown in canonical orientation
    affine = warp_img.affine.copy()
    affine[:3, :3] *= stride
    return np.asarray(nib.as_closest_canonical(nib.Nifti1Image(det, affine)).dataobj)


def _on_fixed_grid(img, fixed, order):
    """Data of img in canonical orientation, resampled to fixed if needed."""
    if img.shape[:3] != fixed.shape[:3] or not np.allclose(img.affine, fixed.affine):
        img = resample_from_to(img, fixed, order=order)
    return np.asarray(img.dataobj)


def _display_window(data, mask):
    """Window and gamma for an image, from its nonzero voxels within mask.

    Stains can have a few very bright regions (e.g. YoPro in the olfactory bulb
    and cerebellum) with most of the brain near the bottom of the range, so the
    gamma puts the median at mid-gray, which keeps the bright regions in range.
    """
    values = data[(data > 0) & mask] if mask is not None else data[data > 0]
    if values.size == 0:
        values = data[data > 0]
    if values.size == 0:
        return 0.0, 1.0, 1.0
    low, median, high = np.percentile(values, [1, 50, 99.5])
    high = max(high, low + 1e-6)
    scaled_median = (median - low) / (high - low)
    gamma = np.log(MEDIAN_GRAY) / np.log(scaled_median) if scaled_median > 0 else 1.0
    return float(low), float(high), float(np.clip(gamma, *GAMMA_RANGE))


def _gray_uri(plane, window):
    low, high, gamma = window
    scaled = np.clip((plane - low) / (high - low), 0, 1) ** gamma
    return _uri(Image.fromarray((scaled * 255).astype(np.uint8)), "JPEG")


def _uri(image, fmt):
    buffer = BytesIO()
    if fmt == "JPEG":
        image.save(buffer, format="JPEG", quality=82)
    else:
        image.save(buffer, format="PNG", optimize=True)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/{fmt.lower()};base64,{encoded}"


def _outline_uri(plane):
    """Transparent PNG with the boundaries between labels."""
    scale = max(1, OUTLINE_WIDTH_PX // max(plane.shape))
    plane = plane.repeat(scale, axis=0).repeat(scale, axis=1)
    edges = np.zeros(plane.shape, dtype=bool)
    edges[:-1, :] |= plane[:-1, :] != plane[1:, :]
    edges[:, :-1] |= plane[:, :-1] != plane[:, 1:]
    rgba = np.zeros(plane.shape + (4,), dtype=np.uint8)
    rgba[edges] = (255, 214, 10, 230)
    return _uri(Image.fromarray(rgba, "RGBA"), "PNG")


def _jacobian_uri(plane):
    scaled = np.log2(np.maximum(plane, 1e-3)) / (2 * LOG_JACOBIAN_RANGE) + 0.5
    rgb = colormaps["RdBu_r"](np.clip(scaled, 0, 1))[..., :3]
    rgb[plane <= 0] = (0.0, 0.0, 0.0)
    return _uri(Image.fromarray((rgb * 255).astype(np.uint8)), "JPEG")


def _to_view(data, axis, index):
    """Slice of a RAS+ array, with superior or anterior up and right on the right."""
    plane = np.take(data, index, axis=axis)
    return np.flipud(plane.T)


def _ncc(a, b, mask):
    a, b = a[mask].astype(np.float64), b[mask].astype(np.float64)
    if a.size < 2 or a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def _jacobian_stats(det, mask):
    values = det[mask] if mask is not None else det.ravel()
    return {
        "median": float(np.median(values)),
        "p1": float(np.percentile(values, 1)),
        "p99": float(np.percentile(values, 99)),
        "folded": float((values <= 0).mean() * 100),
    }


def _slice_indices(mask, axis):
    other = tuple(i for i in range(3) if i != axis)
    present = np.flatnonzero(mask.any(axis=other))
    if present.size == 0:
        present = np.arange(mask.shape[axis])
    return np.linspace(present[0], present[-1], N_SLICES + 2)[1:-1].round().astype(int)


def build_viewer_data(fixed_path, candidates, stains=None):
    """Render the slices for every plane, stain and candidate, and compute metrics."""
    if stains is None:
        stains = [
            {
                "name": "",
                "fixed": fixed_path,
                "warped": [cand["warped"] for cand in candidates],
            }
        ]
    fixed_img = nib.as_closest_canonical(nib.load(stains[0]["fixed"]))
    first_fixed = np.asarray(fixed_img.dataobj, dtype=np.float32)
    zooms = fixed_img.header.get_zooms()[:3]

    labels, dets, stats, masks = [], [], [], []
    for cand in candidates:
        seg = None
        if cand.get("labels"):
            seg = _on_fixed_grid(
                nib.as_closest_canonical(nib.load(cand["labels"])), fixed_img, 0
            ).astype(np.int32)
        labels.append(seg)
        mask = seg > 0 if seg is not None else first_fixed > 0
        masks.append(mask)
        det, det_stats = None, None
        if cand.get("warp"):
            warp_img = nib.load(cand["warp"])
            full = jacobian_determinant(warp_img)
            # the Jacobian layer and mask follow the warp's own grid
            if full.shape == tuple(
                s // JACOBIAN_STRIDE + (s % JACOBIAN_STRIDE > 0)
                for s in first_fixed.shape
            ):
                det = full
                det_stats = _jacobian_stats(
                    det, mask[::JACOBIAN_STRIDE, ::JACOBIAN_STRIDE, ::JACOBIAN_STRIDE]
                )
            else:
                det_stats = _jacobian_stats(full, None)
        dets.append(det)
        stats.append(
            {
                "name": cand["name"],
                "note": cand.get("note", ""),
                "ncc": [],
                "jacobian": det_stats,
            }
        )

    first_labels = next((seg for seg in labels if seg is not None), None)
    if first_labels is not None:
        brain = first_labels > 0
    else:
        # without labels, background would set the window: use the bright half
        brain = (
            first_fixed > np.median(first_fixed[first_fixed > 0])
            if first_fixed.any()
            else first_fixed > 0
        )

    planes = []
    for plane_name, axis, direction in PLANES:
        indices = _slice_indices(brain, axis)
        origin = fixed_img.affine[axis, 3]
        spacing = fixed_img.affine[axis, axis]
        view_shape = _to_view(first_fixed, axis, 0).shape
        width_mm = view_shape[1] * [z for i, z in enumerate(zooms) if i != axis][0]
        height_mm = view_shape[0] * [z for i, z in enumerate(zooms) if i != axis][1]
        plane = {
            "name": plane_name,
            "axis": axis,
            "indices": indices,
            "aspect": f"{width_mm:.4f} / {height_mm:.4f}",
            "positions": [
                f"{direction} {origin + i * spacing:+.2f} mm" for i in indices
            ],
            "fixed": [],
            "candidates": [],
        }
        for seg, det in zip(labels, dets):
            entry = {"warped": [], "warpedFor": [], "labels": None, "jacobian": None}
            if seg is not None:
                entry["labels"] = [
                    _outline_uri(_to_view(seg, axis, i)) for i in indices
                ]
            if det is not None:
                entry["jacobian"] = [
                    _jacobian_uri(_to_view(det, axis, i // JACOBIAN_STRIDE))
                    for i in indices
                ]
            plane["candidates"].append(entry)
        planes.append(plane)

    # the moving image is often the same for every stain, so it is rendered only once
    last_moved = [(None, None)] * len(candidates)
    for stain in stains:
        fixed = _on_fixed_grid(
            nib.as_closest_canonical(nib.load(stain["fixed"])), fixed_img, 1
        ).astype(np.float32)
        fixed_window = _display_window(fixed, brain)
        for plane in planes:
            plane["fixed"].append(
                [
                    _gray_uri(_to_view(fixed, plane["axis"], i), fixed_window)
                    for i in plane["indices"]
                ]
            )
        for c, path in enumerate(stain["warped"]):
            if path is None:
                for plane in planes:
                    plane["candidates"][c]["warpedFor"].append(None)
                stats[c]["ncc"].append(None)
                continue
            is_new = last_moved[c][0] != path
            if is_new:
                moved = _on_fixed_grid(
                    nib.as_closest_canonical(nib.load(path)), fixed_img, 1
                ).astype(np.float32)
                last_moved[c] = (path, moved)
                seg = labels[c]
                moved_window = _display_window(
                    moved, seg > 0 if seg is not None else brain
                )
            moved = last_moved[c][1]
            for plane in planes:
                entry = plane["candidates"][c]
                if is_new:
                    entry["warped"].append(
                        [
                            _gray_uri(_to_view(moved, plane["axis"], i), moved_window)
                            for i in plane["indices"]
                        ]
                    )
                entry["warpedFor"].append(len(entry["warped"]) - 1)
            stats[c]["ncc"].append(_ncc(fixed, moved, masks[c]))

    for plane in planes:
        del plane["axis"], plane["indices"]
    return {
        "stains": [stain["name"] for stain in stains],
        "planes": planes,
        "stats": stats,
    }


def _format(value, spec):
    return "–" if value is None else format(value, spec)


def _stats_table(stats, stain_names):
    rows = []
    for item in stats:
        jac = item["jacobian"] or {}
        folded = jac.get("folded")
        flag = ""
        if folded is not None:
            state = "bad" if folded >= 0.1 else ("warn" if folded > 0 else "good")
            flag = f'<span class="pill {state}">{folded:.3f}%</span>'
        jac_range = f"{jac['p1']:.2f} – {jac['p99']:.2f}" if "p1" in jac else "–"
        note = (
            f'<span class="note">{html.escape(item["note"])}</span>'
            if item["note"]
            else ""
        )
        rows.append(
            f'<tr><th scope="row">{html.escape(item["name"])}{note}</th>'
            + "".join(
                f'<td class="num">{_format(ncc, ".3f")}</td>' for ncc in item["ncc"]
            )
            + f'<td class="num">{_format(jac.get("median"), ".3f")}</td>'
            f'<td class="num">{jac_range}</td>'
            f'<td class="num">{flag or "–"}</td></tr>'
        )
    ncc_headers = "".join(
        f'<th scope="col" class="num">NCC{" " + html.escape(name) if name else " to fixed"}</th>'
        for name in stain_names
    )
    return (
        '<div class="table-wrap"><table><thead><tr><th scope="col">Registration</th>'
        + ncc_headers
        + '<th scope="col" class="num">Median Jacobian</th>'
        '<th scope="col" class="num">Jacobian p1 – p99</th>'
        '<th scope="col" class="num">Folded voxels</th>'
        f"</tr></thead><tbody>{''.join(rows)}</tbody></table></div>"
    )


STYLE = """
:root {
  /* layout: facts and metrics first, then one slice viewer across the page */
  --bg: #f4f6f7; --surface: #ffffff; --ink: #17222b; --muted: #5a6b78;
  --line: #d5dde2; --accent: #0b6e8f; --stage: #0c1216;
  --good: #1a7f4e; --warn: #a66a00; --bad: #b3261e;
  --font-display: "Avenir Next", "Segoe UI Variable Display", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
  --font-body: "Avenir Next", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
  --font-mono: ui-monospace, "SF Mono", "Cascadia Code", Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --bg: #10171c; --surface: #182128; --ink: #e6edf1; --muted: #9db0bc;
    --line: #2d3a43; --accent: #5cc3e4; --stage: #06090b;
    --good: #5fd39a; --warn: #e8b04a; --bad: #ff8a80;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --bg: #10171c; --surface: #182128; --ink: #e6edf1; --muted: #9db0bc;
  --line: #2d3a43; --accent: #5cc3e4; --stage: #06090b;
  --good: #5fd39a; --warn: #e8b04a; --bad: #ff8a80;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
  font: 15px/1.5 var(--font-body); }
main { max-width: 1500px; margin: 0 auto; padding: 28px 20px 48px;
  display: grid; gap: 22px; }
.eyebrow { margin: 0; color: var(--accent); font-weight: 600; font-size: 13px; }
h1 { margin: 2px 0 0; font: 600 27px/1.2 var(--font-display); text-wrap: balance; }
h2 { margin: 0 0 10px; font: 600 17px/1.3 var(--font-display); }
.facts { display: flex; flex-wrap: wrap; gap: 6px 26px; margin: 12px 0 0; }
.facts div { display: flex; gap: 7px; align-items: baseline; }
.facts dt { color: var(--muted); font-size: 13px; }
.facts dd { margin: 0; font-family: var(--font-mono); font-size: 13px;
  overflow-wrap: anywhere; }
section { background: var(--surface); border: 1px solid var(--line);
  border-radius: 8px; padding: 16px 18px; }
.table-wrap { overflow-x: auto; }
table { width: 100%; border-collapse: collapse; font-size: 14px; }
th, td { padding: 7px 12px 7px 0; border-bottom: 1px solid var(--line);
  text-align: left; vertical-align: baseline; }
thead th { color: var(--muted); font-weight: 500; font-size: 13px; }
tbody th { font-weight: 600; }
tbody tr:last-child > * { border-bottom: 0; }
.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.note { display: block; font: 400 12px/1.4 var(--font-mono); color: var(--muted);
  overflow-wrap: anywhere; }
.pill { display: inline-block; padding: 0 8px; border-radius: 10px;
  border: 1px solid currentColor; font-size: 13px; }
.pill.good { color: var(--good); } .pill.warn { color: var(--warn); }
.pill.bad { color: var(--bad); }
.hint { margin: 10px 0 0; color: var(--muted); font-size: 13px; max-width: 80ch; }
.controls { display: flex; flex-wrap: wrap; gap: 12px 22px; align-items: center;
  margin-bottom: 14px; }
.group { display: flex; gap: 2px; padding: 2px; border: 1px solid var(--line);
  border-radius: 7px; }
.group[hidden] { display: none; }
.group button { border: 0; background: transparent; color: var(--ink);
  padding: 5px 11px; border-radius: 5px; font: inherit; font-size: 14px; cursor: pointer; }
.group button:hover { background: color-mix(in srgb, var(--accent) 12%, transparent); }
.group button[aria-pressed="true"] { background: var(--accent); color: var(--surface); }
button:focus-visible, input:focus-visible { outline: 2px solid var(--accent);
  outline-offset: 2px; }
label.check { display: flex; gap: 6px; align-items: center; font-size: 14px; }
.slider { display: flex; gap: 10px; align-items: center; font-size: 14px; flex: 1 1 260px; }
.slider input { flex: 1; accent-color: var(--accent); min-width: 120px; }
.slider output { font-family: var(--font-mono); font-size: 13px; color: var(--muted);
  min-width: 11ch; text-align: right; }
.grid { display: grid; gap: 14px;
  grid-template-columns: repeat(auto-fill, minmax(min(100%, var(--panel, 460px)), 1fr)); }
figure { margin: 0; display: grid; gap: 6px; }
figure[hidden] { display: none; }
figcaption { font-weight: 600; font-size: 14px; }
.stage { position: relative; background: var(--stage); border-radius: 5px;
  overflow: hidden; }
.stage img { display: block; width: 100%; height: 100%; image-rendering: auto; }
.stage img + img { position: absolute; inset: 0; }
.stage img[hidden] { display: none; }
.stage .tag { position: absolute; left: 8px; top: 6px; padding: 1px 7px; border-radius: 4px;
  background: rgb(0 0 0 / 0.62); color: #fff; font-size: 12px; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 18px; margin: 12px 0 0;
  color: var(--muted); font-size: 13px; }
.bar { display: inline-block; width: 84px; height: 9px; border-radius: 2px;
  vertical-align: middle; margin: 0 6px;
  background: linear-gradient(90deg, #053061, #f7f7f7, #67001f); }
.swatch { display: inline-block; width: 18px; height: 3px; vertical-align: middle;
  margin-right: 6px; background: #ffd60a; }
"""

SCRIPT = """
const data = JSON.parse(document.getElementById('viewer-data').textContent);
const state = { plane: 0, stain: 0, slice: 3, layer: data.defaultLayer, labels: true, flicker: null };
const grid = document.getElementById('grid');
const slider = document.getElementById('slice');
const output = document.getElementById('position');
const layerNames = { fixed: 'Fixed', warped: 'Warped moving', jacobian: 'Jacobian' };
const panels = data.stats.map((item) => {
  const figure = document.createElement('figure');
  const caption = document.createElement('figcaption');
  caption.textContent = item.name;
  const stage = document.createElement('div');
  stage.className = 'stage';
  const base = new Image();
  base.alt = '';
  const overlay = new Image();
  overlay.alt = '';
  const tag = document.createElement('span');
  tag.className = 'tag';
  stage.append(base, overlay, tag);
  figure.append(caption, stage);
  grid.append(figure);
  return { figure, stage, base, overlay, tag };
});
function shownLayer() {
  return state.flicker || state.layer;
}
function render() {
  const plane = data.planes[state.plane];
  const i = state.slice;
  output.textContent = plane.positions[i];
  panels.forEach((panel, c) => {
    const cand = plane.candidates[c];
    const warped = cand.warpedFor[state.stain];
    panel.figure.hidden = warped === null;
    if (warped === null) return;
    let layer = shownLayer();
    if (layer === 'jacobian' && !cand.jacobian) layer = 'warped';
    const source = layer === 'fixed' ? plane.fixed[state.stain]
      : layer === 'warped' ? cand.warped[warped] : cand.jacobian;
    panel.stage.style.aspectRatio = plane.aspect;
    panel.base.src = source[i];
    panel.base.alt = `${layerNames[layer]}, ${plane.name}, ${plane.positions[i]}`;
    panel.tag.textContent = layerNames[layer];
    const showLabels = state.labels && cand.labels;
    panel.overlay.hidden = !showLabels;
    if (showLabels) panel.overlay.src = cand.labels[i];
  });
}
function press(group, value) {
  group.querySelectorAll('button').forEach((b) =>
    b.setAttribute('aria-pressed', String(b.dataset.value === String(value))));
}
const stainGroup = document.getElementById('stains');
data.stains.forEach((name, i) => {
  const b = document.createElement('button');
  b.type = 'button';
  b.textContent = name;
  b.dataset.value = i;
  b.addEventListener('click', () => { state.stain = i; press(stainGroup, i); render(); });
  stainGroup.append(b);
});
stainGroup.hidden = data.stains.length < 2;
press(stainGroup, 0);
const planeGroup = document.getElementById('planes');
data.planes.forEach((p, i) => {
  const b = document.createElement('button');
  b.type = 'button';
  b.textContent = p.name;
  b.dataset.value = i;
  b.addEventListener('click', () => { state.plane = i; press(planeGroup, i); render(); });
  planeGroup.append(b);
});
press(planeGroup, 0);
const layerGroup = document.getElementById('layers');
layerGroup.querySelectorAll('button').forEach((b) => b.addEventListener('click', () => {
  state.layer = b.dataset.value;
  press(layerGroup, state.layer);
  render();
}));
press(layerGroup, state.layer);
slider.addEventListener('input', () => { state.slice = Number(slider.value); render(); });
document.getElementById('labels').addEventListener('change', (e) => {
  state.labels = e.target.checked;
  render();
});
let timer = null;
const flicker = document.getElementById('flicker');
flicker.addEventListener('change', () => {
  clearInterval(timer);
  state.flicker = null;
  if (flicker.checked) {
    timer = setInterval(() => {
      state.flicker = state.flicker === 'fixed' ? 'warped' : 'fixed';
      render();
    }, 700);
  }
  render();
});
document.addEventListener('keydown', (e) => {
  if (e.target.tagName === 'INPUT' && e.target.type !== 'range') return;
  if (e.key === 'f') { flicker.click(); }
});
render();
"""


def render_report(
    output_html,
    title,
    eyebrow,
    facts,
    fixed_path,
    candidates,
    hint="",
    panel_width=460,
    default_layer="warped",
    stains=None,
):
    """Write the HTML report for candidates compared against fixed_path.

    default_layer is the layer shown first: fixed when the atlas labels were warped
    along with the moving image, since they then match it by construction.
    stains (see the module docstring) replaces fixed_path and the candidates' warped
    images, to choose between stains in the report.
    """
    viewer = build_viewer_data(fixed_path, candidates, stains)
    viewer["defaultLayer"] = default_layer
    has_labels = any(c.get("labels") for c in candidates)
    has_jacobian = any(entry["jacobian"] for entry in viewer["planes"][0]["candidates"])
    fact_items = "".join(
        f"<div><dt>{html.escape(str(k))}</dt><dd>{html.escape(str(v))}</dd></div>"
        for k, v in facts.items()
    )
    jacobian_button = (
        '<button type="button" data-value="jacobian">Jacobian</button>'
        if has_jacobian
        else ""
    )
    labels_control = (
        '<label class="check"><input type="checkbox" id="labels" checked> '
        "Atlas outlines</label>"
        if has_labels
        else '<input type="checkbox" id="labels" hidden>'
    )
    legend = []
    if has_labels:
        legend.append('<span><span class="swatch"></span>Atlas label boundaries</span>')
    if has_jacobian:
        legend.append(
            f"<span>Jacobian 2<sup>−{LOG_JACOBIAN_RANGE:g}</sup> shrink"
            f'<span class="bar"></span>2<sup>{LOG_JACOBIAN_RANGE:g}</sup> expand; '
            "black is folded</span>"
        )
    payload = json.dumps(viewer, separators=(",", ":")).replace("</", "<\\/")
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>{STYLE}</style>
</head>
<body>
<main style="--panel: {panel_width}px">
<header>
<p class="eyebrow">{html.escape(eyebrow)}</p>
<h1>{html.escape(title)}</h1>
<dl class="facts">{fact_items}</dl>
</header>
<section aria-labelledby="metrics-title">
<h2 id="metrics-title">Metrics</h2>
{_stats_table(viewer["stats"], viewer["stains"])}
<p class="hint">{html.escape(hint)} These numbers flag problems but do not rank registrations: judge alignment from the slices. NCC is computed within the atlas or the fixed image's nonzero voxels, and is only a rough guide between different contrasts. Jacobians are within the same region, and include any affine scaling. Each image has its own display window and gamma, so brightness is not comparable between images.</p>
</section>
<section aria-labelledby="viewer-title">
<h2 id="viewer-title">Slices</h2>
<div class="controls">
<div class="group" id="stains" role="group" aria-label="Stain"></div>
<div class="group" id="planes" role="group" aria-label="Plane"></div>
<div class="group" id="layers" role="group" aria-label="Layer">
<button type="button" data-value="fixed">Fixed</button>
<button type="button" data-value="warped">Warped moving</button>
{jacobian_button}
</div>
{labels_control}
<label class="check"><input type="checkbox" id="flicker"> Flicker fixed and warped (F)</label>
<label class="slider" for="slice">Slice
<input type="range" id="slice" min="0" max="{N_SLICES - 1}" value="3" step="1">
<output id="position" for="slice"></output></label>
</div>
<div class="grid" id="grid"></div>
<div class="legend">{"".join(legend)}</div>
</section>
</main>
<script type="application/json" id="viewer-data">{payload}</script>
<script>{SCRIPT}</script>
</body>
</html>
"""
    with open(output_html, "w") as f:
        f.write(page)
