"""Otsu threshold sweep QC report.

For a given stain and multiotsu segmentation method, sweeps over a range of
threshold values (spanning the 1st–99th percentile of the image) and generates
2D crops at multiple positions, producing a self-contained HTML report.  The
report can be visually assessed to select the optimal threshold before running
the full segmentation pipeline.

The otsu histogram PNG figures produced by the ``multiotsu`` rule are embedded
in the report alongside the sweep visualizations.

This is a Snakemake script; the ``snakemake`` object is automatically provided
when executed as part of a Snakemake workflow.
"""

import base64
from io import BytesIO

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


def main():
    zarrnii_kwargs = snakemake.params.zarrnii_kwargs
    n_thresholds = snakemake.params.n_thresholds
    n_crops = snakemake.params.n_crops
    patch_size = snakemake.params.patch_size
    level = snakemake.params.level
    pct_lo, pct_hi = snakemake.params.hist_percentile_range

    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc
    subject = snakemake.wildcards.subject

    # ------------------------------------------------------------------
    # Load a downsampled pyramid level for speed
    # ------------------------------------------------------------------
    print("Loading image for threshold sweep QC...")
    znimg = None
    for ds_offset in [3, 2, 1, 0]:
        try:
            candidate = ZarrNii.from_ome_zarr(
                snakemake.input.corrected, level=level + ds_offset, **zarrnii_kwargs
            )
            znimg = candidate
            print(f"  Loaded at pyramid level {level + ds_offset}")
            break
        except Exception as exc:
            print(f"  Level {level + ds_offset} not available: {exc}")

    if znimg is None:
        znimg = ZarrNii.from_ome_zarr(snakemake.input.corrected, **zarrnii_kwargs)

    # ------------------------------------------------------------------
    # Get array data — shape is (C, Z, Y, X) or (Z, Y, X)
    # ------------------------------------------------------------------
    data = znimg.data.compute()
    if data.ndim == 4:
        data = data[0]  # drop channel dim → (Z, Y, X)

    # ------------------------------------------------------------------
    # Compute percentile-based range for the sweep
    # ------------------------------------------------------------------
    flat = data.ravel().astype(np.float32)
    range_lo = float(np.percentile(flat, pct_lo))
    range_hi = float(np.percentile(flat, pct_hi))
    print(
        f"  Sweep range [{pct_lo}th–{pct_hi}th percentile]: "
        f"[{range_lo:.3f}, {range_hi:.3f}]"
    )

    # ------------------------------------------------------------------
    # Generate evenly-spaced threshold values across the percentile range
    # ------------------------------------------------------------------
    thresholds = np.linspace(range_lo, range_hi, n_thresholds)

    # ------------------------------------------------------------------
    # Select crop positions: evenly spaced along Y at mid-Z
    # ------------------------------------------------------------------
    Z, Y, X = data.shape
    z_mid = Z // 2
    half_patch = patch_size // 2

    crop_boxes = []
    for i in range(n_crops):
        y_pos = int(Y * (i + 0.5) / n_crops)
        x_pos = X // 2
        y0 = max(0, y_pos - half_patch)
        y1 = min(Y, y_pos + half_patch)
        x0 = max(0, x_pos - half_patch)
        x1 = min(X, x_pos + half_patch)
        crop_boxes.append((z_mid, y0, y1, x0, x1))

    # Pre-normalise the raw image crops (shared across all thresholds)
    img_crops_norm = [
        _norm(data[z, y0:y1, x0:x1], range_lo, range_hi)
        for (z, y0, y1, x0, x1) in crop_boxes
    ]

    # ------------------------------------------------------------------
    # Build a figure for each threshold value
    # ------------------------------------------------------------------
    threshold_entries = []
    for thresh in thresholds:
        fig, axes = plt.subplots(1, n_crops, figsize=(n_crops * 3, 3))
        fig.suptitle(f"Threshold = {thresh:.1f}", fontsize=10, fontweight="bold")
        if n_crops == 1:
            axes = [axes]

        for ax, (z, y0, y1, x0, x1), crop_norm in zip(axes, crop_boxes, img_crops_norm):
            mask_crop = (data[z, y0:y1, x0:x1] > thresh).astype(np.float32)
            ax.imshow(crop_norm, cmap="gray")
            mask_ma = np.ma.masked_where(mask_crop < 0.5, mask_crop)
            ax.imshow(mask_ma, cmap="Reds", alpha=0.6, vmin=0, vmax=1)
            ax.set_xticks([])
            ax.set_yticks([])

        plt.tight_layout()
        threshold_entries.append(
            {"threshold": float(thresh), "b64": _fig_to_base64(fig)}
        )

    # ------------------------------------------------------------------
    # Embed the otsu histogram PNG from the multiotsu rule
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

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Otsu Threshold Sweep QC &ndash; {subject}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
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
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #3498db;
            padding-left: 10px;
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
        <h1>Otsu Threshold Sweep QC Report</h1>

        <div class="info-box">
            <p><strong>Subject:</strong> {subject}</p>
            <p><strong>Stain:</strong> {stain}</p>
            <p><strong>Method:</strong> {desc}</p>
            <p><strong>Intensity range ({pct_lo}th&ndash;{pct_hi}th percentile):</strong>
               {range_lo:.1f} &ndash; {range_hi:.1f}</p>
            <p><strong>Number of threshold values:</strong> {n_thresholds}</p>
        </div>

        <h2>Otsu Histogram (Multi-level Otsu Thresholds)</h2>
        <p>The histogram below shows the multi-Otsu threshold positions computed
        from the bias-field corrected image.  Use these as a guide when selecting
        the optimal threshold from the sweep below.</p>
        <div class="figure">
            <img src="{otsu_hist_b64}" alt="Otsu Histogram">
        </div>

        <h2>Threshold Sweep (2D Crops)</h2>
        <p>Each row shows {n_crops} axial crops at different Y positions (mid-Z
        slice) with the resulting binary mask overlaid in red.  Use this to select
        a threshold that captures the signal of interest without excessive
        background segmentation.</p>

        {sweep_html}

    </div>
</body>
</html>"""

    with open(snakemake.output.html, "w") as fh:
        fh.write(html)

    print(f"Saved threshold sweep QC to {snakemake.output.html}")


if __name__ == "__main__":
    main()
