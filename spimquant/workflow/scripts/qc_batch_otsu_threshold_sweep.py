"""Batch Otsu threshold sweep QC report.

For a given stain and multiotsu segmentation method, uses the batch-wide
Multi-Otsu thresholds (computed from the aggregated histogram across all
subjects) to show one mid-volume 2D crop per subject at each threshold
value.  This enables visual assessment of whether the same threshold works
across an entire acquisition batch.

The aggregated histogram PNG is embedded at the top of the report, followed
by a sweep of threshold values, each showing one patch per subject side by
side.

This is a Snakemake script; the ``snakemake`` object is automatically provided
when executed as part of a Snakemake workflow.
"""

import base64
import re
from io import BytesIO
from pathlib import Path

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from zarrnii import ZarrNii


def _fig_to_base64(fig):
    """Convert a matplotlib figure to a base64-encoded PNG data URI."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return f"data:image/png;base64,{b64}"


def _file_to_base64_png(path):
    """Read a PNG file and return a base64 data URI."""
    with open(path, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode("utf-8")
    return f"data:image/png;base64,{b64}"


def _norm(arr, lo, hi):
    """Linearly normalise *arr* to [0, 1] using the given bounds."""
    if hi > lo:
        return np.clip((arr.astype(np.float32) - lo) / (hi - lo), 0.0, 1.0)
    return np.zeros_like(arr, dtype=np.float32)


def _subject_label(zarr_path):
    """Extract a short subject label from a file path."""
    m = re.search(r"sub-([^_/\\]+)", zarr_path)
    if m:
        return f"sub-{m.group(1)}"
    return Path(zarr_path).stem


def main():
    zarrnii_kwargs = snakemake.params.zarrnii_kwargs
    n_thresholds = snakemake.params.n_thresholds
    patch_size = snakemake.params.patch_size
    level = snakemake.params.level
    pct_lo, pct_hi = snakemake.params.hist_percentile_range

    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc

    zarrs = snakemake.input.zarrs
    n_subjects = len(zarrs)
    subject_labels = [_subject_label(p) for p in zarrs]

    print(f"Batch threshold sweep QC: stain={stain}, desc={desc}")
    print(f"  {n_subjects} subject(s): {subject_labels}")

    # ------------------------------------------------------------------
    # Load one mid-volume crop per subject
    # ------------------------------------------------------------------
    half_patch = patch_size // 2
    subject_crops = []
    subject_ranges = []

    for zarr_path, label in zip(zarrs, subject_labels):
        print(f"  Loading crop from {label} …")
        znimg = None
        for ds_offset in [5,4,3, 2, 1, 0]:
            try:
                candidate = ZarrNii.from_ome_zarr(
                    zarr_path, level=level + ds_offset, **zarrnii_kwargs
                )
                znimg = candidate
                print(f"    Loaded at pyramid level {level + ds_offset}")
                break
            except Exception as exc:
                print(f"    Level {level + ds_offset} not available: {exc}")

        if znimg is None:
            znimg = ZarrNii.from_ome_zarr(zarr_path, **zarrnii_kwargs)

        data = znimg.data.compute()
        if data.ndim == 4:
            data = data[0]  # drop channel dim → (Z, Y, X)

        flat = data.ravel().astype(np.float32)
        range_lo = float(np.percentile(flat, pct_lo))
        range_hi = float(np.percentile(flat, pct_hi))
        subject_ranges.append((range_lo, range_hi))

        Z, Y, X = data.shape
        z_mid = Z // 2
        y_pos = Y // 2
        x_pos = X // 2
        y0 = max(0, y_pos - half_patch)
        y1 = min(Y, y_pos + half_patch)
        x0 = max(0, x_pos - half_patch)
        x1 = min(X, x_pos + half_patch)
        subject_crops.append(data[z_mid, y0:y1, x0:x1])

    # ------------------------------------------------------------------
    # Determine global sweep range and threshold values
    # ------------------------------------------------------------------
    global_lo = float(np.min([r[0] for r in subject_ranges]))
    global_hi = float(np.max([r[1] for r in subject_ranges]))
    thresholds = np.linspace(global_lo, global_hi, n_thresholds)
    print(f"  Sweep range: [{global_lo:.3f}, {global_hi:.3f}]")

    # ------------------------------------------------------------------
    # Build one figure per threshold: columns = subjects
    # ------------------------------------------------------------------
    threshold_entries = []
    for thresh in thresholds:
        fig, axes = plt.subplots(
            1, n_subjects, figsize=(max(3 * n_subjects, 6), 3), squeeze=False
        )
        axes = axes[0]
        fig.suptitle(f"Threshold = {thresh:.1f}", fontsize=10, fontweight="bold")

        for ax, crop, (rlo, rhi), label in zip(
            axes, subject_crops, subject_ranges, subject_labels
        ):
            crop_norm = _norm(crop, rlo, rhi)
            mask_crop = (crop > thresh).astype(np.float32)
            ax.imshow(crop_norm, cmap="gray")
            mask_ma = np.ma.masked_where(mask_crop < 0.5, mask_crop)
            ax.imshow(mask_ma, cmap="Reds", alpha=0.6, vmin=0, vmax=1)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(label, fontsize=7)

        plt.tight_layout()
        threshold_entries.append(
            {"threshold": float(thresh), "b64": _fig_to_base64(fig)}
        )

    # ------------------------------------------------------------------
    # Embed the batch Otsu histogram PNG
    # ------------------------------------------------------------------
    otsu_hist_b64 = _file_to_base64_png(snakemake.input.thresholds_png)

    # ------------------------------------------------------------------
    # Generate HTML report
    # ------------------------------------------------------------------
    sweep_html_parts = []
    for entry in threshold_entries:
        thresh_label = f"{entry['threshold']:.1f}"
        sweep_html_parts.append(
            f'<div class="figure">'
            f"<h3>Threshold = {thresh_label}</h3>"
            f'<img src="{entry["b64"]}" alt="threshold {thresh_label}">'
            f"</div>"
        )
    sweep_html = "\n".join(sweep_html_parts)

    subjects_str = ", ".join(subject_labels)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Batch Otsu Threshold Sweep QC &ndash; {stain} / {desc}</title>
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
            border-bottom: 3px solid #e74c3c;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #e74c3c;
            padding-left: 10px;
        }}
        .info-box {{
            background-color: #ecf0f1;
            border-left: 4px solid #e74c3c;
            padding: 15px;
            margin: 20px 0;
        }}
        .info-box p {{
            margin: 5px 0;
        }}
        .figure {{
            margin: 20px 0;
            text-align: center;
        }}
        .figure img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Batch Otsu Threshold Sweep QC Report</h1>

        <div class="info-box">
            <p><strong>Stain:</strong> {stain}</p>
            <p><strong>Method:</strong> {desc}</p>
            <p><strong>Subjects ({n_subjects}):</strong> {subjects_str}</p>
            <p><strong>Intensity sweep range
               ({pct_lo}th&ndash;{pct_hi}th percentile):</strong>
               {global_lo:.1f} &ndash; {global_hi:.1f}</p>
            <p><strong>Number of threshold values:</strong> {n_thresholds}</p>
        </div>

        <h2>Batch Histogram (Multi-level Otsu Thresholds from Aggregated Data)</h2>
        <p>The histogram below was computed from the <strong>aggregated</strong>
        intensity data across all subjects in this batch.  The vertical lines
        show the Multi-Otsu threshold positions.  These batch-wide thresholds
        are used as a fixed reference for all subjects so that segmentation is
        consistent across acquisitions.</p>
        <div class="figure">
            <img src="{otsu_hist_b64}" alt="Batch Otsu Histogram">
        </div>

        <h2>Threshold Sweep (one crop per subject)</h2>
        <p>Each row corresponds to one threshold value.  Within each row, every
        column shows a mid-volume 2D crop from one subject with the resulting
        binary mask overlaid in red.  Use this to identify a threshold that
        reliably captures the signal of interest across <em>all</em> subjects
        without excessive background segmentation.</p>

        {sweep_html}

    </div>
</body>
</html>"""

    with open(snakemake.output.html, "w") as fh:
        fh.write(html)

    print(f"Saved batch threshold sweep QC to {snakemake.output.html}")


if __name__ == "__main__":
    main()
