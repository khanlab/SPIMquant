"""Dataset QC report: one self-contained HTML page on all scans.

Takes the per-scan outputs of qc_scan.py. Scans may have different stains: the
report has a column per stain (labels that differ only in case, or are listed
as aliases, share one), and a scan without a stain has an empty cell there.

Metrics are flagged when they are far from the other scans of the same stain,
which needs a few scans to mean anything, or beyond a fixed limit (saturation,
field-of-view cut-off, atlas overlap). Flags point at what to look at; they do
not decide whether a scan is usable.
"""

import html
import json

import numpy as np
import pandas as pd

# key, label, what a problem looks like ("low", "high" or "both"), help
STAIN_METRICS = [
    ("snr", "SNR", "low", "Tissue median above background, over background noise."),
    (
        "tissue_to_background",
        "Tissue / background",
        "low",
        "Tissue median over background median.",
    ),
    (
        "dynamic_range",
        "Signal range",
        "both",
        "Brightest tissue (99.9th percentile) over the tissue median: high for "
        "sparse bright labelling, near 1 for flat signal.",
    ),
    ("p50", "Tissue median", "both", "Median intensity in tissue."),
    ("p99_9", "Tissue p99.9", "both", "99.9th percentile of intensity in tissue."),
    (
        "profile_cv",
        "Depth non-uniformity",
        "high",
        "Variation of the median tissue intensity from slice to slice, along the "
        "worst axis: high with attenuation or uneven clearing.",
    ),
    (
        "sharpness",
        "Sharpness",
        "low",
        "Fine detail relative to coarse structure in the crops (median); "
        "compare between scans of the same stain.",
    ),
    (
        "saturated",
        "Saturated",
        "high",
        "Fraction of pixels at the top of the data type, in the worst crop.",
    ),
]
SCAN_METRICS = [
    ("tissue_volume_mm3", "Tissue volume (mm³)", "both", "Volume of the tissue mask."),
    (
        "edge_touch",
        "At field-of-view edge",
        "high",
        "Fraction of the tissue surface on the faces of the image: high when the "
        "sample is cut off.",
    ),
    (
        "reg_dice",
        "Atlas overlap",
        "low",
        "Dice overlap of the registered atlas with the tissue mask.",
    ),
    (
        "reg_ncc",
        "Template similarity",
        "low",
        "Correlation of the registered scan with the template.",
    ),
]
# fixed limits: metric -> (direction, limit)
LIMITS = {
    "saturated": ("high", 0.001),
    "edge_touch": ("high", 0.05),
    "reg_dice": ("low", 0.8),
}
OUTLIER_Z = 3.5
MIN_SCANS_FOR_OUTLIERS = 4
ENTITY_TAGS = {
    "subject": "sub",
    "session": "ses",
    "sample": "sample",
    "acquisition": "acq",
}


def canonical_stains(labels, aliases=None):
    """Map stain labels to the name of their column in the report.

    Labels listed as aliases take their canonical name; labels that differ only
    in case share the first spelling seen.
    """
    alias_of = {}
    for canonical, names in (aliases or {}).items():
        alias_of[canonical.lower()] = canonical
        for name in names:
            alias_of[name.lower()] = canonical
    mapping, by_lower = {}, {}
    for label in labels:
        key = label.lower()
        mapping[label] = alias_of.get(key) or by_lower.setdefault(key, label)
    return mapping


def scan_label(entities, varying):
    """BIDS-style name of a scan from the entities that tell the scans apart."""
    return "_".join(
        f"{tag}-{entities[name]}"
        for name, tag in ENTITY_TAGS.items()
        if name in entities and (name == "subject" or name in varying)
    )


def flag(value, others, direction, key):
    """Why a value is flagged, or None.

    Against the fixed limit of the metric if it has one, and against the other
    scans (robust z-score beyond OUTLIER_Z in the direction that is a problem).
    """
    if value is None:
        return None
    if key in LIMITS:
        limit_direction, limit = LIMITS[key]
        if (value > limit) if limit_direction == "high" else (value < limit):
            return f"{'above' if limit_direction == 'high' else 'below'} {limit:g}"
    values = np.array([v for v in others if v is not None], dtype=float)
    if values.size < MIN_SCANS_FOR_OUTLIERS:
        return None
    median = np.median(values)
    spread = 1.4826 * np.median(np.abs(values - median))
    if spread == 0:
        return None
    z = (value - median) / spread
    if (direction in ("high", "both") and z > OUTLIER_Z) or (
        direction in ("low", "both") and z < -OUTLIER_Z
    ):
        return f"{'high' if z > 0 else 'low'} for this stain (z = {z:+.1f})"
    return None


def build_report_data(
    scans, title, stain_aliases=None, participants=None, group_by=None
):
    """The data of the report, from the per-scan outputs of qc_scan.py."""
    all_labels = [s for scan in scans for s in scan["stains"]]
    canonical = canonical_stains(all_labels, stain_aliases)
    stains = list(dict.fromkeys(canonical[s] for s in all_labels))

    varying = {
        name
        for name in ENTITY_TAGS
        if len({scan["entities"].get(name) for scan in scans}) > 1
    }
    group_columns = []
    by_participant = {}
    if participants is not None and "participant_id" in participants.columns:
        group_columns = [
            c
            for c in participants.columns
            if c != "participant_id" and (not group_by or c in group_by)
        ]
        by_participant = participants.set_index("participant_id").to_dict("index")

    out_scans = []
    for scan in scans:
        entities = scan["entities"]
        registration = scan["registration"] or {}
        row = by_participant.get(f"sub-{entities.get('subject')}", {})
        out_scans.append(
            {
                "id": scan_label(entities, set(ENTITY_TAGS)),
                "label": scan_label(entities, varying),
                "groups": {c: str(row.get(c, "")) for c in group_columns},
                "registration": scan["registration"],
                "meta": scan["meta"],
                "planes": scan["planes"],
                "regions": scan["regions"],
                "correlation": scan["scan"]["correlation"],
                "stainNames": [canonical[s] for s in scan["stains"]],
                "metrics": {
                    "tissue_volume_mm3": scan["scan"]["tissue_volume_mm3"],
                    "edge_touch": scan["scan"]["edge_touch"],
                    "reg_dice": registration.get("dice"),
                    "reg_ncc": registration.get("ncc"),
                },
                "flags": {},
                "channels": {
                    canonical[s]: {**scan["channels"][s], "label": s, "flags": {}}
                    for s in scan["stains"]
                },
            }
        )
    out_scans.sort(key=lambda s: s["id"])

    for key, _, direction, _ in SCAN_METRICS:
        values = [s["metrics"][key] for s in out_scans]
        for scan in out_scans:
            reason = flag(scan["metrics"][key], values, direction, key)
            if reason:
                scan["flags"][key] = reason
    windows, max_hi = {}, {}
    for stain in stains:
        having = [s["channels"][stain] for s in out_scans if stain in s["channels"]]
        for key, _, direction, _ in STAIN_METRICS:
            values = [c["metrics"].get(key) for c in having]
            for channel in having:
                reason = flag(channel["metrics"].get(key), values, direction, key)
                if reason:
                    channel["flags"][key] = reason
        # one display window for all scans of the stain
        windows[stain] = [
            float(np.median([c["window"][0] for c in having])),
            float(np.median([c["window"][1] for c in having])),
        ]
        max_hi[stain] = float(max(c["hi"] for c in having))

    region_names = list(
        dict.fromkeys(r["name"] for s in out_scans for r in s["regions"])
    )
    return {
        "title": title,
        "stains": stains,
        "groupColumns": group_columns,
        "regionNames": region_names,
        "stainMetrics": [
            {"key": k, "label": l, "help": h} for k, l, _, h in STAIN_METRICS
        ],
        "scanMetrics": [
            {"key": k, "label": l, "help": h} for k, l, _, h in SCAN_METRICS
        ],
        "windows": windows,
        "maxHi": max_hi,
        "scans": out_scans,
    }


def metrics_table(data):
    """Long table of the metrics: one row per scan and stain."""
    rows = []
    for scan in data["scans"]:
        base = {"scan": scan["id"], **scan["groups"]}
        base.update(scan["metrics"])
        base["scan_flags"] = "; ".join(f"{k}: {v}" for k, v in scan["flags"].items())
        for stain, channel in scan["channels"].items():
            row = {**base, "stain": stain, "stain_label": channel["label"]}
            for metric in data["stainMetrics"]:
                row[metric["key"]] = channel["metrics"].get(metric["key"])
            row["stain_flags"] = "; ".join(
                f"{k}: {v}" for k, v in channel["flags"].items()
            )
            rows.append(row)
    return pd.DataFrame(rows)


STYLE = """
:root {
  --bg: #f4f6f7; --surface: #ffffff; --ink: #17222b; --muted: #5a6b78;
  --line: #d5dde2; --accent: #0b6e8f; --stage: #0c1216;
  --good: #1a7f4e; --warn: #a66a00; --bad: #b3261e;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a;
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
    --s1: #3987e5; --s2: #d95926; --s3: #199e70;
  }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--ink);
  font: 15px/1.5 var(--font-body); }
main { max-width: 1600px; margin: 0 auto; padding: 28px 20px 48px;
  display: grid; gap: 22px; }
.eyebrow { margin: 0; color: var(--accent); font-weight: 600; font-size: 13px; }
h1 { margin: 2px 0 0; font: 600 27px/1.2 var(--font-display); text-wrap: balance; }
h2 { margin: 0 0 10px; font: 600 17px/1.3 var(--font-display); }
h3 { margin: 18px 0 8px; font: 600 15px/1.3 var(--font-display); }
.facts { display: flex; flex-wrap: wrap; gap: 6px 26px; margin: 12px 0 0; }
.facts div { display: flex; gap: 7px; align-items: baseline; }
.facts dt { color: var(--muted); font-size: 13px; }
.facts dd { margin: 0; font-family: var(--font-mono); font-size: 13px;
  overflow-wrap: anywhere; }
section { background: var(--surface); border: 1px solid var(--line);
  border-radius: 8px; padding: 16px 18px; min-width: 0; }
.table-wrap { overflow-x: auto; }
table { border-collapse: collapse; font-size: 14px; }
th, td { padding: 7px 12px 7px 0; border-bottom: 1px solid var(--line);
  text-align: left; vertical-align: top; }
thead th { color: var(--muted); font-weight: 500; font-size: 13px; }
tbody th { font-weight: 600; }
.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.sub { display: block; font: 400 12px/1.4 var(--font-mono); color: var(--muted); }
.flag { color: var(--warn); font-size: 12px; display: block; }
.flag::before { content: "\\26A0\\FE0E  "; }
.hint { margin: 10px 0 0; color: var(--muted); font-size: 13px; max-width: 90ch; }
.controls { display: flex; flex-wrap: wrap; gap: 10px 20px; align-items: center;
  margin-bottom: 14px; font-size: 14px; }
.controls label { display: flex; gap: 6px; align-items: center; }
select, input[type="text"] { font: inherit; font-size: 14px; color: var(--ink);
  background: var(--surface); border: 1px solid var(--line); border-radius: 6px;
  padding: 4px 8px; }
.group { display: flex; flex-wrap: wrap; gap: 2px; padding: 2px;
  border: 1px solid var(--line); border-radius: 7px; }
.group button, button.plain { border: 0; background: transparent; color: var(--ink);
  padding: 5px 11px; border-radius: 5px; font: inherit; font-size: 14px; cursor: pointer; }
button.plain { border: 1px solid var(--line); }
.group button:hover, button.plain:hover {
  background: color-mix(in srgb, var(--accent) 12%, transparent); }
.group button[aria-pressed="true"] { background: var(--accent); color: var(--surface); }
button:focus-visible, input:focus-visible, select:focus-visible, td.cell:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 2px; }
input[type="range"] { accent-color: var(--accent); width: 150px; }
.window { position: sticky; top: 0; z-index: 5; background: var(--surface);
  border: 1px solid var(--line); border-radius: 8px; padding: 10px 18px;
  display: flex; flex-wrap: wrap; gap: 8px 20px; align-items: center; font-size: 14px; }
.window label { display: flex; gap: 6px; align-items: center; }
.window output { font-family: var(--font-mono); font-size: 13px; color: var(--muted);
  min-width: 6ch; }
.stage { position: relative; background: var(--stage); border-radius: 5px;
  overflow: hidden; line-height: 0; }
.stage canvas { display: block; width: 100%; height: auto; }
.stage img { position: absolute; inset: 0; width: 100%; height: 100%; }
.stage img[hidden] { display: none; }
#matrix tbody th { min-width: 150px; }
#matrix td.cell { cursor: pointer; width: 170px; min-width: 130px; padding-right: 14px; }
#matrix td.cell:hover .stage { outline: 2px solid var(--accent); }
#matrix td.missing { color: var(--muted); background: repeating-linear-gradient(
  135deg, transparent 0 7px, color-mix(in srgb, var(--line) 55%, transparent) 7px 8px); }
#matrix .value { font-variant-numeric: tabular-nums; font-size: 13px; margin-top: 4px;
  display: block; }
.grid { display: grid; gap: 14px;
  grid-template-columns: repeat(auto-fill, minmax(min(100%, var(--panel, 300px)), 1fr)); }
figure { margin: 0; display: grid; gap: 5px; align-content: start; min-width: 0; }
figcaption { font-size: 13px; font-weight: 600; overflow-wrap: anywhere; }
figcaption .sub { font-weight: 400; }
.crops { display: grid; gap: 10px; align-items: start; justify-content: start; }
.crops figure { max-width: 512px; }
.crops .head { font-size: 13px; color: var(--muted); font-weight: 500; }
.crops .region { font-size: 13px; font-weight: 600; align-self: center; }
.two { display: grid; gap: 22px; grid-template-columns: repeat(auto-fit, minmax(min(100%, 420px), 1fr)); }
svg.chart { width: 100%; height: auto; display: block; }
svg.chart text { font: 11px var(--font-body); fill: var(--muted); }
.legend { display: flex; flex-wrap: wrap; gap: 4px 16px; font-size: 13px;
  color: var(--muted); margin-top: 6px; }
.key { display: inline-block; width: 14px; height: 3px; vertical-align: middle;
  margin-right: 6px; border-radius: 2px; }
.rating { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; margin-top: 12px; }
.rating input[type="text"] { flex: 1 1 240px; }
.pill { display: inline-block; padding: 0 7px; border-radius: 10px; font-size: 12px;
  border: 1px solid currentColor; margin-top: 3px; }
.pill.pass { color: var(--good); } .pill.check { color: var(--warn); }
.pill.fail { color: var(--bad); }
#tip { position: fixed; pointer-events: none; z-index: 20; background: var(--ink);
  color: var(--surface); padding: 3px 8px; border-radius: 4px; font-size: 12px;
  white-space: nowrap; }
#tip[hidden] { display: none; }
"""

SCRIPT = r"""
const D = JSON.parse(document.getElementById('qc-data').textContent);
const $ = (id) => document.getElementById(id);
const el = (tag, props = {}, ...kids) => {
  const node = Object.assign(document.createElement(tag), props);
  node.append(...kids.filter((k) => k !== null && k !== undefined));
  return node;
};
const fmt = (v) => v === null || v === undefined ? '–'
  : Math.abs(v) >= 1000 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1)
  : Math.abs(v) >= 0.01 || v === 0 ? v.toFixed(2) : v.toExponential(1);
const state = {
  metric: D.stainMetrics[0].key, plane: 0, scan: 0, stain: D.stains[0], slice: 0,
  outlines: true, shared: true, gamma: 0.6, compareStain: D.stains[0], compareWhat: '',
  compareGroup: '', windows: JSON.parse(JSON.stringify(D.windows)),
};
const slices = D.scans[0].planes[0].outlines.length;
state.slice = Math.floor(slices / 2);

/* ---- images: stored 8-bit with a square-root transfer of [0, hi]; decoded to
   intensities and shown with the window of their stain ---- */
const decoded = new Map();
function decode(uri) {
  if (!decoded.has(uri)) {
    decoded.set(uri, new Promise((resolve) => {
      const img = new Image();
      img.onload = () => {
        const c = document.createElement('canvas');
        c.width = img.naturalWidth; c.height = img.naturalHeight;
        const ctx = c.getContext('2d');
        ctx.drawImage(img, 0, 0);
        const rgba = ctx.getImageData(0, 0, c.width, c.height).data;
        const gray = new Uint8Array(c.width * c.height);
        for (let i = 0; i < gray.length; i++) gray[i] = rgba[4 * i];
        resolve({ w: c.width, h: c.height, gray });
      };
      img.src = uri;
    }));
  }
  return decoded.get(uri);
}
function windowOf(stain, channel) {
  return state.shared ? state.windows[stain] : channel.window;
}
async function paint(canvas) {
  const { uri, stain, channel } = canvas.spec;
  const im = await decode(uri);
  if (canvas.spec.uri !== uri) return;
  const [lo, hi] = windowOf(stain, channel);
  const lut = new Uint8Array(256);
  for (let v = 0; v < 256; v++) {
    const x = channel.hi * (v / 255) ** 2;
    const t = Math.min(1, Math.max(0, (x - lo) / Math.max(hi - lo, 1e-9)));
    lut[v] = Math.round(255 * t ** state.gamma);
  }
  if (canvas.width !== im.w) canvas.width = im.w;
  if (canvas.height !== im.h) canvas.height = im.h;
  const ctx = canvas.getContext('2d');
  const out = ctx.createImageData(im.w, im.h);
  for (let i = 0; i < im.gray.length; i++) {
    const g = lut[im.gray[i]];
    out.data[4 * i] = g; out.data[4 * i + 1] = g; out.data[4 * i + 2] = g;
    out.data[4 * i + 3] = 255;
  }
  ctx.putImageData(out, 0, 0);
}
function stage(uri, stain, channel, outline, alt) {
  const canvas = el('canvas', { role: 'img', ariaLabel: alt });
  canvas.spec = { uri, stain, channel };
  const box = el('div', { className: 'stage' }, canvas);
  if (outline) {
    const img = el('img', { src: outline, alt: '', className: 'outline' });
    img.hidden = !state.outlines;
    box.append(img);
  }
  paint(canvas);
  return box;
}
let pending = false;
function repaint() {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => {
    pending = false;
    document.querySelectorAll('canvas').forEach((c) => { if (c.spec) paint(c); });
  });
}
function buttons(group, items, current, onPick) {
  group.replaceChildren(...items.map(([value, label]) => {
    const b = el('button', { type: 'button', textContent: label });
    b.setAttribute('aria-pressed', String(value === current));
    b.addEventListener('click', () => onPick(value));
    return b;
  }));
}
function options(select, items, current) {
  select.replaceChildren(...items.map(([value, label]) =>
    el('option', { value, textContent: label, selected: value === current })));
}

/* ---- ratings, kept in this browser ---- */
const storeKey = 'spimquant-qc:' + D.title;
let ratings = {};
try { ratings = JSON.parse(localStorage.getItem(storeKey)) || {}; } catch (e) {}
function saveRatings() {
  try { localStorage.setItem(storeKey, JSON.stringify(ratings)); } catch (e) {}
}
const ratingNames = { pass: 'Pass', check: 'Check', fail: 'Fail' };

/* ---- window controls ---- */
const toSlider = (x, stain) => Math.round(1000 * Math.sqrt(Math.max(0, x) / D.maxHi[stain]));
const fromSlider = (s, stain) => D.maxHi[stain] * (s / 1000) ** 2;
function syncWindow() {
  const [lo, hi] = state.windows[state.stain];
  $('w-stain').textContent = state.stain;
  $('w-lo').value = toSlider(lo, state.stain); $('w-hi').value = toSlider(hi, state.stain);
  $('w-lo-out').textContent = fmt(lo); $('w-hi-out').textContent = fmt(hi);
  $('w-lo').disabled = $('w-hi').disabled = !state.shared;
}
['w-lo', 'w-hi'].forEach((id, i) => $(id).addEventListener('input', (e) => {
  const w = state.windows[state.stain];
  w[i] = fromSlider(Number(e.target.value), state.stain);
  if (w[1] <= w[0]) w[1 - i] = i ? w[1] * 0.9 : w[0] * 1.1 + 1e-6;
  syncWindow(); repaint(); drawHistogram();
}));
$('w-gamma').addEventListener('input', (e) => {
  state.gamma = Number(e.target.value); $('w-gamma-out').textContent = state.gamma.toFixed(2);
  repaint();
});
$('w-shared').addEventListener('change', (e) => {
  state.shared = e.target.checked; syncWindow(); repaint(); drawHistogram();
});
$('w-reset').addEventListener('click', () => {
  state.windows[state.stain] = [...D.windows[state.stain]];
  syncWindow(); repaint(); drawHistogram();
});
function setStain(stain) {
  state.stain = stain; state.compareStain = stain;
  syncWindow();
}

/* ---- matrix: scans by stains ---- */
function renderMatrix() {
  const head = el('tr', {}, el('th', { scope: 'col', textContent: 'Scan' }),
    ...D.stains.map((s) => el('th', { scope: 'col', textContent: s })));
  const rows = D.scans.map((scan, i) => {
    const th = el('th', { scope: 'row' }, scan.label);
    const groups = Object.values(scan.groups).filter(Boolean).join(' · ');
    if (groups) th.append(el('span', { className: 'sub', textContent: groups }));
    if (!scan.registration) th.append(el('span', { className: 'sub', textContent: 'not registered' }));
    Object.entries(scan.flags).forEach(([key, why]) => th.append(el('span', {
      className: 'flag',
      textContent: D.scanMetrics.find((m) => m.key === key).label + ' ' + why })));
    const cells = D.stains.map((stain) => {
      const ch = scan.channels[stain];
      if (!ch) return el('td', { className: 'missing', textContent: 'not stained' });
      const td = el('td', { className: 'cell', tabIndex: 0 });
      td.append(stage(ch.overview[state.plane][state.slice], stain, ch, null,
        `${scan.label} ${stain}`));
      const name = D.stainMetrics.find((m) => m.key === state.metric).label;
      td.append(el('span', { className: 'value',
        textContent: `${name} ${fmt(ch.metrics[state.metric])}` }));
      if (ch.label !== stain) td.append(el('span', { className: 'sub', textContent: ch.label }));
      Object.entries(ch.flags).forEach(([key, why]) => td.append(el('span', {
        className: 'flag',
        textContent: D.stainMetrics.find((m) => m.key === key).label + ' ' + why })));
      const r = ratings[scan.id + '|' + stain];
      if (r && r.rating) td.append(el('span', { className: 'pill ' + r.rating,
        textContent: ratingNames[r.rating] }));
      const open = () => {
        state.scan = i; setStain(stain); renderDetail(); renderCompare();
        $('detail').scrollIntoView({ behavior: 'smooth' });
      };
      td.addEventListener('click', open);
      td.addEventListener('keydown', (e) => { if (e.key === 'Enter') open(); });
      return td;
    });
    return el('tr', {}, th, ...cells);
  });
  $('matrix').replaceChildren(el('thead', {}, head), el('tbody', {}, ...rows));
}

/* ---- one scan ---- */
function renderDetail() {
  const scan = D.scans[state.scan];
  if (!scan.channels[state.stain]) setStain(scan.stainNames[0]);
  const stain = state.stain, ch = scan.channels[stain];
  options($('d-scan'), D.scans.map((s, i) => [String(i), s.label]), String(state.scan));
  buttons($('d-stains'), scan.stainNames.map((s) => [s, s]), stain, (s) => {
    setStain(s); renderDetail(); renderCompare();
  });
  const meta = scan.meta, reg = scan.registration;
  const facts = {
    'Size (x, y, z)': meta.shape.join(' × '),
    'Voxel (µm)': meta.voxel_um.join(' × '),
    'Orientation': meta.orientation,
    'Bits': meta.bits,
    'Registered': reg ? `${reg.stain} to ${reg.template}, atlas ${reg.atlas}` : 'no',
    ...scan.groups,
  };
  $('d-facts').replaceChildren(...Object.entries(facts).map(([k, v]) =>
    el('div', {}, el('dt', { textContent: k }), el('dd', { textContent: String(v) }))));

  $('d-slice').max = slices - 1; $('d-slice').value = state.slice;
  $('d-outlines').closest('label').hidden = !reg;
  $('d-overview').replaceChildren(...scan.planes.map((plane, p) => el('figure', {},
    el('figcaption', { textContent: plane.name }),
    stage(ch.overview[p][state.slice], stain, ch, plane.outlines[state.slice],
      `${scan.label} ${stain} ${plane.name}`))));

  const grid = $('d-crops');
  grid.style.gridTemplateColumns =
    `minmax(70px, 110px) repeat(${scan.stainNames.length}, minmax(0, 512px))`;
  const cells = [el('span'), ...scan.stainNames.map((s) =>
    el('span', { className: 'head', textContent: s }))];
  scan.regions.forEach((region, r) => {
    cells.push(el('span', { className: 'region', textContent: region.name }));
    scan.stainNames.forEach((s) => {
      const c = scan.channels[s], crop = c.crops[r];
      if (!crop) { cells.push(el('span', { className: 'sub', textContent: 'outside the scan' })); return; }
      const sat = crop.saturated ? ` · saturated ${(100 * crop.saturated).toFixed(2)}%` : '';
      cells.push(el('figure', {}, stage(crop.image, s, c, null, `${scan.label} ${s} ${region.name}`),
        el('span', { className: 'sub', textContent: `sharpness ${fmt(crop.sharpness)}${sat}` })));
    });
  });
  grid.replaceChildren(...cells);
  $('d-crops-hint').textContent = `One plane at pyramid level ${meta.crop_level}, at the same place in every stain` +
    (reg ? ', in the largest or chosen atlas regions.' : '; spread through the tissue, as this scan is not registered.');

  const head = el('tr', {}, el('th', { scope: 'col', textContent: 'Metric' }),
    ...scan.stainNames.map((s) => el('th', { scope: 'col', className: 'num', textContent: s })));
  const rows = D.stainMetrics.map((m) => el('tr', {},
    el('th', { scope: 'row', textContent: m.label, title: m.help }),
    ...scan.stainNames.map((s) => {
      const c = scan.channels[s];
      const td = el('td', { className: 'num', textContent: fmt(c.metrics[m.key]) });
      if (c.flags[m.key]) td.append(el('span', { className: 'flag', textContent: c.flags[m.key] }));
      return td;
    })));
  $('d-metrics').replaceChildren(el('thead', {}, head), el('tbody', {}, ...rows));
  $('d-scanmetrics').replaceChildren(el('tbody', {}, ...D.scanMetrics.map((m) => {
    const td = el('td', { className: 'num', textContent: fmt(scan.metrics[m.key]) });
    if (scan.flags[m.key]) td.append(el('span', { className: 'flag', textContent: scan.flags[m.key] }));
    return el('tr', {}, el('th', { scope: 'row', textContent: m.label, title: m.help }), td);
  })));
  const corr = scan.correlation;
  $('d-corr').replaceChildren(...(corr && scan.stainNames.length > 1 ? [
    el('thead', {}, el('tr', {}, el('th'), ...scan.stainNames.map((s) =>
      el('th', { scope: 'col', className: 'num', textContent: s })))),
    el('tbody', {}, ...scan.stainNames.map((s, i) => el('tr', {},
      el('th', { scope: 'row', textContent: s }),
      ...corr[i].map((v, j) => el('td', { className: 'num', textContent: i === j ? '' : fmt(v) })))))] : []));
  $('d-corr-title').hidden = !(corr && scan.stainNames.length > 1);

  drawHistogram(); drawProfiles();
  const key = scan.id + '|' + stain, r = ratings[key] || {};
  $('r-title').textContent = `Rating of ${scan.label}, ${stain}`;
  buttons($('r-buttons'), Object.entries(ratingNames), r.rating, (value) => {
    ratings[key] = { ...ratings[key], rating: r.rating === value ? '' : value };
    saveRatings(); renderDetail(); renderMatrix();
  });
  $('r-note').value = r.note || '';
  $('r-note').onchange = (e) => {
    ratings[key] = { ...ratings[key], note: e.target.value }; saveRatings();
  };
}

/* ---- charts ---- */
const SVG = 'http://www.w3.org/2000/svg';
const svg = (tag, attrs, text) => {
  const node = document.createElementNS(SVG, tag);
  Object.entries(attrs).forEach(([k, v]) => node.setAttribute(k, v));
  if (text !== undefined) node.textContent = text;
  return node;
};
const tip = $('tip');
function showTip(e, text) {
  tip.textContent = text; tip.hidden = false;
  tip.style.left = Math.min(e.clientX + 12, innerWidth - tip.offsetWidth - 8) + 'px';
  tip.style.top = (e.clientY - 30) + 'px';
}
function drawHistogram() {
  const scan = D.scans[state.scan], ch = scan.channels[state.stain], h = ch.histogram;
  const chart = $('d-hist');
  chart.replaceChildren();
  if (!h) return;
  const W = 420, H = 150, L = 8, R = 8, T = 8, B = 26;
  const n = h.counts.length, peak = Math.max(...h.counts), bw = (W - L - R) / n;
  const x = (v) => L + (W - L - R) * Math.log(Math.max(v, h.lo) / h.lo) / Math.log(h.hi / h.lo);
  const edge = (i) => h.lo * (h.hi / h.lo) ** (i / n);
  h.counts.forEach((c, i) => {
    const bh = (H - T - B) * c / peak;
    const bar = svg('rect', { x: L + i * bw + 0.5, y: H - B - bh, width: Math.max(bw - 1, 1),
      height: bh, fill: 'var(--s1)' });
    const text = `${fmt(edge(i))} – ${fmt(edge(i + 1))}: ${c.toLocaleString()} voxels`;
    const hit = svg('rect', { x: L + i * bw, y: T, width: bw, height: H - T - B, fill: 'transparent' });
    hit.addEventListener('pointermove', (e) => showTip(e, text));
    hit.addEventListener('pointerleave', () => { tip.hidden = true; });
    chart.append(bar, hit);
  });
  chart.append(svg('line', { x1: L, x2: W - R, y1: H - B, y2: H - B, stroke: 'var(--line)' }));
  windowOf(state.stain, ch).forEach((v, i) => {
    if (v < h.lo || v > h.hi) return;
    chart.append(svg('line', { x1: x(v), x2: x(v), y1: T, y2: H - B, stroke: 'var(--ink)',
      'stroke-dasharray': '3 3' }));
    chart.append(svg('text', { x: x(v), y: H - B + 12, 'text-anchor': i ? 'end' : 'start' },
      (i ? 'white ' : 'black ') + fmt(v)));
  });
  chart.append(svg('text', { x: L, y: H - 2 }, fmt(h.lo)));
  chart.append(svg('text', { x: W - R, y: H - 2, 'text-anchor': 'end' }, fmt(h.hi)));
  chart.append(svg('text', { x: W / 2, y: H - 2, 'text-anchor': 'middle' }, 'intensity in tissue (log scale)'));
}
function drawProfiles() {
  const scan = D.scans[state.scan], ch = scan.channels[state.stain];
  const chart = $('d-prof');
  chart.replaceChildren();
  const W = 420, H = 150, L = 30, R = 8, T = 8, B = 22;
  const all = ch.profiles.flat().filter((v) => v !== null);
  if (!all.length) return;
  const top = Math.max(1.2, Math.max(...all) * 1.05);
  const y = (v) => T + (H - T - B) * (1 - v / top);
  [0, 1].forEach((v) => {
    chart.append(svg('line', { x1: L, x2: W - R, y1: y(v), y2: y(v), stroke: 'var(--line)' }));
    chart.append(svg('text', { x: L - 5, y: y(v) + 4, 'text-anchor': 'end' }, v));
  });
  chart.append(svg('text', { x: L, y: H - 4 }, 'first slice'));
  chart.append(svg('text', { x: W - R, y: H - 4, 'text-anchor': 'end' }, 'last slice'));
  ch.profiles.forEach((profile, a) => {
    const px = (i) => L + (W - L - R) * i / Math.max(profile.length - 1, 1);
    let d = '', pen = false;
    profile.forEach((v, i) => {
      if (v === null) { pen = false; return; }
      d += `${pen ? 'L' : 'M'}${px(i).toFixed(1)} ${y(v).toFixed(1)}`; pen = true;
    });
    chart.append(svg('path', { d, fill: 'none', stroke: `var(--s${a + 1})`, 'stroke-width': 2,
      'stroke-linejoin': 'round' }));
  });
  const hit = svg('rect', { x: L, y: T, width: W - L - R, height: H - T - B, fill: 'transparent' });
  hit.addEventListener('pointermove', (e) => {
    const box = hit.getBoundingClientRect(), f = (e.clientX - box.left) / box.width;
    showTip(e, ch.profiles.map((p, a) => {
      const v = p[Math.round(f * (p.length - 1))];
      return `${scan.planes[a].name} ${v === null ? '–' : v.toFixed(2)}`;
    }).join(' · '));
  });
  hit.addEventListener('pointerleave', () => { tip.hidden = true; });
  chart.append(hit);
  $('d-prof-legend').replaceChildren(...scan.planes.map((p, a) => el('span', {},
    el('span', { className: 'key', style: `background: var(--s${a + 1})` }),
    `through ${p.name} slices`)));
}

/* ---- one stain across scans ---- */
function renderCompare() {
  buttons($('c-stains'), D.stains.map((s) => [s, s]), state.compareStain, (s) => {
    setStain(s); renderCompare(); if (D.scans[state.scan].channels[s]) renderDetail();
  });
  options($('c-what'), [['', 'Overview'], ...D.regionNames.map((n) => [n, 'Crop: ' + n])],
    state.compareWhat);
  options($('c-group'), [['', 'Scan'], ...D.groupColumns.map((c) => [c, c])], state.compareGroup);
  $('c-group').closest('label').hidden = !D.groupColumns.length;
  const stain = state.compareStain, g = state.compareGroup;
  const having = D.scans.map((scan, i) => ({ scan, i })).filter(({ scan }) => scan.channels[stain]);
  if (g) having.sort((a, b) => a.scan.groups[g].localeCompare(b.scan.groups[g]) || a.i - b.i);
  $('c-grid').replaceChildren(...having.map(({ scan, i }) => {
    const ch = scan.channels[stain];
    let box = null, note = '';
    if (state.compareWhat) {
      const r = scan.regions.findIndex((x) => x.name === state.compareWhat);
      const crop = r < 0 ? null : ch.crops[r];
      if (crop) {
        box = stage(crop.image, stain, ch, null, `${scan.label} ${stain} ${state.compareWhat}`);
        note = `sharpness ${fmt(crop.sharpness)}`;
      }
    } else {
      box = stage(ch.overview[state.plane][state.slice], stain, ch,
        scan.planes[state.plane].outlines[state.slice], `${scan.label} ${stain}`);
      note = `SNR ${fmt(ch.metrics.snr)} · signal range ${fmt(ch.metrics.dynamic_range)}`;
    }
    const caption = el('figcaption', {}, (g ? scan.groups[g] + ' · ' : '') + scan.label,
      el('span', { className: 'sub', textContent: box ? note : 'no crop of this region' }));
    const figure = el('figure', {}, caption, box);
    figure.style.cursor = 'pointer';
    figure.addEventListener('click', () => {
      state.scan = i; setStain(stain); renderDetail();
      $('detail').scrollIntoView({ behavior: 'smooth' });
    });
    return figure;
  }));
  $('c-hint').textContent = `${having.length} of ${D.scans.length} scans have ${stain}.`;
}

/* ---- wiring ---- */
options($('m-metric'), D.stainMetrics.map((m) => [m.key, m.label]), state.metric);
$('m-metric').addEventListener('change', (e) => { state.metric = e.target.value; renderMatrix(); });
function renderPlanes() {
  const items = D.scans[0].planes.map((p, i) => [i, p.name]);
  const pick = (i) => { state.plane = i; renderPlanes(); renderMatrix(); renderCompare(); };
  buttons($('m-planes'), items, state.plane, pick);
  buttons($('c-planes'), items, state.plane, pick);
}
$('d-scan').addEventListener('change', (e) => { state.scan = Number(e.target.value); renderDetail(); });
$('d-slice').addEventListener('input', (e) => {
  state.slice = Number(e.target.value); renderDetail(); renderMatrix(); renderCompare();
});
$('d-outlines').addEventListener('change', (e) => {
  state.outlines = e.target.checked;
  document.querySelectorAll('img.outline').forEach((img) => { img.hidden = !state.outlines; });
});
$('c-what').addEventListener('change', (e) => { state.compareWhat = e.target.value; renderCompare(); });
$('c-group').addEventListener('change', (e) => { state.compareGroup = e.target.value; renderCompare(); });
$('r-export').addEventListener('click', () => {
  const lines = ['scan\tstain\trating\tnote'];
  D.scans.forEach((scan) => Object.keys(scan.channels).forEach((stain) => {
    const r = ratings[scan.id + '|' + stain] || {};
    lines.push([scan.id, stain, r.rating || '', (r.note || '').replace(/[\t\n]/g, ' ')].join('\t'));
  }));
  const a = el('a', { download: 'qcratings.tsv',
    href: URL.createObjectURL(new Blob([lines.join('\n') + '\n'], { type: 'text/tab-separated-values' })) });
  a.click();
});
syncWindow(); renderPlanes(); renderMatrix(); renderDetail(); renderCompare();
"""


def render_report(output_html, data):
    """Write the report as one HTML file with everything embedded."""
    scans = data["scans"]
    registered = sum(1 for s in scans if s["registration"])
    facts = {
        "Scans": len(scans),
        "Stains": ", ".join(data["stains"]),
        "Registered": f"{registered} of {len(scans)}",
    }
    fact_items = "".join(
        f"<div><dt>{html.escape(str(k))}</dt><dd>{html.escape(str(v))}</dd></div>"
        for k, v in facts.items()
    )
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    title = html.escape(data["title"])
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>QC: {title}</title>
<style>{STYLE}</style>
</head>
<body>
<main>
<header>
<p class="eyebrow">SPIMquant dataset QC</p>
<h1>{title}</h1>
<dl class="facts">{fact_items}</dl>
</header>
<div class="window" role="group" aria-label="Display window">
<span>Window of <strong id="w-stain"></strong></span>
<label>Black <input type="range" id="w-lo" min="0" max="1000"><output id="w-lo-out"></output></label>
<label>White <input type="range" id="w-hi" min="0" max="1000"><output id="w-hi-out"></output></label>
<label>Gamma <input type="range" id="w-gamma" min="0.3" max="1" step="0.05" value="0.6"><output id="w-gamma-out">0.60</output></label>
<label><input type="checkbox" id="w-shared" checked> Same window for all scans</label>
<button type="button" class="plain" id="w-reset">Reset</button>
</div>
<section aria-labelledby="matrix-title">
<h2 id="matrix-title">Scans by stain</h2>
<div class="controls">
<label>Metric <select id="m-metric"></select></label>
<div class="group" id="m-planes" role="group" aria-label="Plane"></div>
</div>
<div class="table-wrap"><table id="matrix"></table></div>
<p class="hint">Select a cell to open that scan and stain below. A warning marks a metric that is far from the other scans of the same stain (robust z-score beyond {OUTLIER_Z:g}, with at least {MIN_SCANS_FOR_OUTLIERS} scans) or beyond a fixed limit; it says where to look, not whether the scan is usable. Every scan of a stain is shown with the same window unless that is switched off above, so brightness is comparable down a column.</p>
</section>
<section aria-labelledby="detail-title" id="detail">
<h2 id="detail-title">One scan</h2>
<div class="controls">
<label>Scan <select id="d-scan"></select></label>
<div class="group" id="d-stains" role="group" aria-label="Stain"></div>
<label>Slice <input type="range" id="d-slice" min="0" step="1"></label>
<label><input type="checkbox" id="d-outlines" checked> Atlas outlines</label>
</div>
<dl class="facts" id="d-facts"></dl>
<h3>Overview</h3>
<div class="grid" id="d-overview" style="--panel: 340px"></div>
<h3>Full-resolution crops</h3>
<div class="crops" id="d-crops"></div>
<p class="hint" id="d-crops-hint"></p>
<div class="two">
<div>
<h3>Metrics</h3>
<div class="table-wrap"><table id="d-metrics"></table></div>
<div class="table-wrap"><table id="d-scanmetrics" style="margin-top: 14px"></table></div>
<h3 id="d-corr-title">Correlation between stains</h3>
<div class="table-wrap"><table id="d-corr"></table></div>
</div>
<div>
<h3>Intensity histogram</h3>
<svg class="chart" id="d-hist" viewBox="0 0 420 150" role="img" aria-label="Histogram of the intensity in tissue"></svg>
<h3>Intensity through the sample</h3>
<svg class="chart" id="d-prof" viewBox="0 0 420 150" role="img" aria-label="Median tissue intensity by slice, relative to the tissue median"></svg>
<div class="legend" id="d-prof-legend"></div>
<p class="hint">Median tissue intensity of each slice, relative to the median of all tissue. A slope or dip is attenuation or uneven clearing along that direction.</p>
</div>
</div>
<h3 id="r-title">Rating</h3>
<div class="rating">
<div class="group" id="r-buttons" role="group" aria-label="Rating"></div>
<input type="text" id="r-note" placeholder="Note" aria-label="Note">
<button type="button" class="plain" id="r-export">Download ratings (TSV)</button>
</div>
<p class="hint">Ratings and notes are kept in this browser only; download them to keep or share them.</p>
</section>
<section aria-labelledby="compare-title">
<h2 id="compare-title">One stain across scans</h2>
<div class="controls">
<div class="group" id="c-stains" role="group" aria-label="Stain"></div>
<label>Show <select id="c-what"></select></label>
<div class="group" id="c-planes" role="group" aria-label="Plane"></div>
<label>Order by <select id="c-group"></select></label>
</div>
<div class="grid" id="c-grid" style="--panel: 280px"></div>
<p class="hint" id="c-hint"></p>
</section>
</main>
<div id="tip" hidden></div>
<script type="application/json" id="qc-data">{payload}</script>
<script>{SCRIPT}</script>
</body>
</html>
"""
    with open(output_html, "w") as f:
        f.write(page)
