"""N4 bias field correction parameter tuning – HTML report.

Aggregates the per-combination QC PNGs produced by ``qc_n4_tune_png`` and
assembles them into a self-contained, scrollable HTML report.

Report structure
----------------
1. **Summary comparison table** – one cell per (spline_spacing, iters)
   combination.  Each cell shows key statistics for the corrected image and
   bias field (min, max, mean, p50, p99, bias-field CV) so different parameter
   settings can be compared at a glance.
2. **Detailed rows** – one row per combination, with the embedded PNG (base64)
   and a per-image statistics sub-table (uncorrected, bias field, corrected).

Inputs
------
The ``snakemake`` object exposes:

* ``snakemake.input.pngs``       – list of PNG paths (expand order: spline ×
                                   iters, spline iterates slowest).
* ``snakemake.input.corrected``  – matching list of corrected NIfTI paths.
* ``snakemake.input.biasfield``  – matching list of bias-field NIfTI paths.
* ``snakemake.input.nii``        – single uncorrected NIfTI path.
* ``snakemake.params.spline_spacings``  – list of spline-spacing strings.
* ``snakemake.params.iters_list``       – list of iterations strings.

This is a Snakemake script; the ``snakemake`` object is automatically provided
when executed as part of a Snakemake workflow.
"""

import base64
import itertools

import nibabel as nib
import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _file_to_base64_png(path):
    """Read a PNG file and return a ``data:image/png;base64,…`` URI."""
    with open(path, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def _load_stats(path):
    """Load a NIfTI and return a summary statistics dict."""
    img = nib.load(path)
    data = np.asarray(img.dataobj).ravel().astype(np.float32)
    pctls = np.percentile(data, [1, 5, 25, 50, 75, 95, 99])
    mean = float(data.mean())
    std = float(data.std())
    return {
        "min": float(data.min()),
        "max": float(data.max()),
        "mean": mean,
        "std": std,
        "cv": float(std / mean) if mean != 0.0 else float("nan"),
        "p1": float(pctls[0]),
        "p5": float(pctls[1]),
        "p25": float(pctls[2]),
        "p50": float(pctls[3]),
        "p75": float(pctls[4]),
        "p95": float(pctls[5]),
        "p99": float(pctls[6]),
    }


def _fmt(value, decimals=3):
    """Format a float for display in the HTML table."""
    if isinstance(value, float) and (value != value):  # NaN
        return "N/A"
    return f"{value:.{decimals}f}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    spline_spacings = snakemake.params.spline_spacings
    iters_list = snakemake.params.iters_list
    subject = snakemake.wildcards.subject
    stain = snakemake.wildcards.stain

    combos = list(itertools.product(spline_spacings, iters_list))
    n_combos = len(combos)

    print(
        f"Building N4 tuning report for subject={subject}, stain={stain}, "
        f"{n_combos} parameter combinations…"
    )

    # -----------------------------------------------------------------------
    # Load data for all combinations
    # -----------------------------------------------------------------------
    png_paths = list(snakemake.input.pngs)
    corrected_paths = list(snakemake.input.corrected)
    biasfield_paths = list(snakemake.input.biasfield)
    raw_path = snakemake.input.nii

    assert len(png_paths) == n_combos, (
        f"Expected {n_combos} PNGs, got {len(png_paths)}"
    )
    assert len(corrected_paths) == n_combos
    assert len(biasfield_paths) == n_combos

    print("  Computing statistics for uncorrected image…")
    raw_stats = _load_stats(raw_path)

    combo_data = []
    for idx, (spline, iters) in enumerate(combos):
        print(f"  Loading stats for spline={spline}, iters={iters}…")
        combo_data.append(
            {
                "spline": spline,
                "iters": iters,
                "png_b64": _file_to_base64_png(png_paths[idx]),
                "corrected_stats": _load_stats(corrected_paths[idx]),
                "biasfield_stats": _load_stats(biasfield_paths[idx]),
            }
        )

    # -----------------------------------------------------------------------
    # Build summary comparison table (2-D grid: rows=spline, cols=iters)
    # -----------------------------------------------------------------------
    iters_header = "".join(
        f"<th>Iters={i}</th>" for i in iters_list
    )
    summary_rows_html = []
    for spline in spline_spacings:
        cells = []
        for iters in iters_list:
            entry = next(
                d for d in combo_data if d["spline"] == spline and d["iters"] == iters
            )
            cs = entry["corrected_stats"]
            bf = entry["biasfield_stats"]
            cell = (
                f"<td>"
                f"<strong>Corrected</strong><br>"
                f"min: {_fmt(cs['min'])}<br>"
                f"max: {_fmt(cs['max'])}<br>"
                f"p50: {_fmt(cs['p50'])}<br>"
                f"p99: {_fmt(cs['p99'])}<br>"
                f"mean: {_fmt(cs['mean'])}<br>"
                f"<br><strong>Bias Field</strong><br>"
                f"min: {_fmt(bf['min'])}<br>"
                f"max: {_fmt(bf['max'])}<br>"
                f"mean: {_fmt(bf['mean'])}<br>"
                f"CV: {_fmt(bf['cv'])}<br>"
                f"p50: {_fmt(bf['p50'])}<br>"
                f"p99: {_fmt(bf['p99'])}"
                f"</td>"
            )
            cells.append(cell)
        row_html = (
            f"<tr><th>Spline={spline} mm</th>{''.join(cells)}</tr>"
        )
        summary_rows_html.append(row_html)

    summary_table_html = f"""
<table class="summary-table">
  <thead>
    <tr><th>Parameter</th>{iters_header}</tr>
  </thead>
  <tbody>
    {''.join(summary_rows_html)}
  </tbody>
</table>"""

    # -----------------------------------------------------------------------
    # Build raw-image stats row for reference
    # -----------------------------------------------------------------------
    raw_stats_html = (
        f"<p><strong>Uncorrected image stats:</strong> "
        f"min={_fmt(raw_stats['min'])} | max={_fmt(raw_stats['max'])} | "
        f"mean={_fmt(raw_stats['mean'])} | "
        f"p50={_fmt(raw_stats['p50'])} | p99={_fmt(raw_stats['p99'])}</p>"
    )

    # -----------------------------------------------------------------------
    # Build detailed rows (one per combination)
    # -----------------------------------------------------------------------
    detail_rows_html = []
    for entry in combo_data:
        spline = entry["spline"]
        iters = entry["iters"]
        cs = entry["corrected_stats"]
        bf = entry["biasfield_stats"]

        detail_stat_table = f"""
<table class="stat-table">
  <thead>
    <tr>
      <th>Metric</th>
      <th>Uncorrected</th>
      <th>Bias Field</th>
      <th>Corrected</th>
    </tr>
  </thead>
  <tbody>
    <tr><td>min</td><td>{_fmt(raw_stats['min'])}</td>
        <td>{_fmt(bf['min'])}</td><td>{_fmt(cs['min'])}</td></tr>
    <tr><td>max</td><td>{_fmt(raw_stats['max'])}</td>
        <td>{_fmt(bf['max'])}</td><td>{_fmt(cs['max'])}</td></tr>
    <tr><td>mean</td><td>{_fmt(raw_stats['mean'])}</td>
        <td>{_fmt(bf['mean'])}</td><td>{_fmt(cs['mean'])}</td></tr>
    <tr><td>std</td><td>{_fmt(raw_stats['std'])}</td>
        <td>{_fmt(bf['std'])}</td><td>{_fmt(cs['std'])}</td></tr>
    <tr><td>CV</td><td>{_fmt(raw_stats['cv'])}</td>
        <td>{_fmt(bf['cv'])}</td><td>{_fmt(cs['cv'])}</td></tr>
    <tr><td>p1</td><td>{_fmt(raw_stats['p1'])}</td>
        <td>{_fmt(bf['p1'])}</td><td>{_fmt(cs['p1'])}</td></tr>
    <tr><td>p5</td><td>{_fmt(raw_stats['p5'])}</td>
        <td>{_fmt(bf['p5'])}</td><td>{_fmt(cs['p5'])}</td></tr>
    <tr><td>p25</td><td>{_fmt(raw_stats['p25'])}</td>
        <td>{_fmt(bf['p25'])}</td><td>{_fmt(cs['p25'])}</td></tr>
    <tr><td>p50</td><td>{_fmt(raw_stats['p50'])}</td>
        <td>{_fmt(bf['p50'])}</td><td>{_fmt(cs['p50'])}</td></tr>
    <tr><td>p75</td><td>{_fmt(raw_stats['p75'])}</td>
        <td>{_fmt(bf['p75'])}</td><td>{_fmt(cs['p75'])}</td></tr>
    <tr><td>p95</td><td>{_fmt(raw_stats['p95'])}</td>
        <td>{_fmt(bf['p95'])}</td><td>{_fmt(cs['p95'])}</td></tr>
    <tr><td>p99</td><td>{_fmt(raw_stats['p99'])}</td>
        <td>{_fmt(bf['p99'])}</td><td>{_fmt(cs['p99'])}</td></tr>
  </tbody>
</table>"""

        detail_rows_html.append(
            f"""
<div class="combo-row" id="combo-spline{spline}-iters{iters}">
  <h2>Spline spacing: {spline} mm &mdash; Iterations per level: {iters}</h2>
  <img src="{entry['png_b64']}"
       alt="N4 tuning QC: spline={spline}, iters={iters}">
  <h3>Quantitative statistics</h3>
  {detail_stat_table}
</div>"""
        )

    details_html = "\n".join(detail_rows_html)

    # -----------------------------------------------------------------------
    # Assemble full HTML
    # -----------------------------------------------------------------------
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>N4 Parameter Tuning QC &ndash; {subject} / {stain}</title>
  <style>
    body {{
      font-family: Arial, sans-serif;
      margin: 20px;
      background-color: #f5f5f5;
    }}
    .container {{
      max-width: 1600px;
      margin: 0 auto;
      background-color: white;
      padding: 30px;
      box-shadow: 0 0 10px rgba(0,0,0,0.1);
    }}
    h1 {{
      color: #2c3e50;
      border-bottom: 3px solid #3498db;
      padding-bottom: 10px;
    }}
    h2 {{
      color: #2c3e50;
      margin-top: 40px;
      border-left: 4px solid #3498db;
      padding-left: 10px;
      background-color: #eaf4fb;
      padding: 8px 12px;
    }}
    h3 {{
      color: #34495e;
      margin-top: 20px;
    }}
    .info-box {{
      background-color: #ecf0f1;
      border-left: 4px solid #3498db;
      padding: 15px;
      margin: 20px 0;
    }}
    .info-box p {{
      margin: 5px 0;
    }}
    /* Summary table */
    .summary-table {{
      border-collapse: collapse;
      width: 100%;
      margin: 20px 0;
      font-size: 12px;
    }}
    .summary-table th, .summary-table td {{
      border: 1px solid #ccc;
      padding: 8px 12px;
      text-align: left;
      vertical-align: top;
    }}
    .summary-table thead th {{
      background-color: #2c3e50;
      color: white;
      font-weight: bold;
    }}
    .summary-table tbody tr:nth-child(odd) {{
      background-color: #f9f9f9;
    }}
    .summary-table tbody th {{
      background-color: #ecf0f1;
      font-weight: bold;
    }}
    /* Per-combination stat table */
    .stat-table {{
      border-collapse: collapse;
      width: 100%;
      margin: 10px 0;
      font-size: 12px;
    }}
    .stat-table th, .stat-table td {{
      border: 1px solid #ddd;
      padding: 5px 10px;
      text-align: right;
    }}
    .stat-table thead th {{
      background-color: #3498db;
      color: white;
      text-align: center;
    }}
    .stat-table tbody tr:nth-child(odd) {{
      background-color: #f9f9f9;
    }}
    .stat-table tbody td:first-child {{
      text-align: left;
      font-weight: bold;
    }}
    /* Per-combination visual row */
    .combo-row {{
      margin: 30px 0;
      padding: 20px;
      border: 1px solid #ddd;
      border-radius: 4px;
      background-color: #fafafa;
    }}
    .combo-row img {{
      max-width: 100%;
      height: auto;
      display: block;
      margin: 10px auto;
      border: 1px solid #ccc;
    }}
    /* Jump links */
    .jump-links {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin: 15px 0;
    }}
    .jump-links a {{
      background-color: #3498db;
      color: white;
      padding: 4px 10px;
      border-radius: 3px;
      text-decoration: none;
      font-size: 12px;
    }}
    .jump-links a:hover {{
      background-color: #2980b9;
    }}
  </style>
</head>
<body>
  <div class="container">
    <h1>N4 Bias Field Correction &mdash; Parameter Tuning QC Report</h1>

    <div class="info-box">
      <p><strong>Subject:</strong> {subject}</p>
      <p><strong>Stain:</strong> {stain}</p>
      <p><strong>Spline spacing values (mm):</strong> {", ".join(spline_spacings)}</p>
      <p><strong>Iterations per level:</strong> {", ".join(iters_list)}</p>
      <p><strong>Total combinations:</strong> {n_combos}</p>
    </div>

    {raw_stats_html}

    <h2>Summary Comparison Table</h2>
    <p>Each cell shows statistics for the <em>corrected</em> image and the
    estimated <em>bias field</em> for that parameter combination.  A lower
    bias-field coefficient of variation (CV) indicates a smoother field;
    higher spline spacing produces a smoother (coarser) field.  Use this
    table together with the visual comparisons below to select the best
    parameters.</p>
    {summary_table_html}

    <h2>Quick Navigation</h2>
    <div class="jump-links">
      {''.join(
          f'<a href="#combo-spline{d["spline"]}-iters{d["iters"]}">'
          f'Spline {d["spline"]} mm / Iters {d["iters"]}</a>'
          for d in combo_data
      )}
    </div>

    <h2>Detailed Visual Comparisons</h2>
    <p>Each row below shows the uncorrected image, estimated bias field, and
    corrected image for one parameter combination (axial, coronal, and sagittal
    middle slices, plus intensity histograms and statistics).  Scroll down to
    compare parameter settings.</p>

    {details_html}

  </div>
</body>
</html>"""

    with open(snakemake.output.html, "w") as fh:
        fh.write(html)

    print(f"Saved N4 tuning report to {snakemake.output.html}")


if __name__ == "__main__":
    main()
