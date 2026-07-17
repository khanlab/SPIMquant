"""N4 bias field correction parameter tuning – per-combination QC figure.

Generates a multi-panel PNG for a single (spline_spacing, iters) parameter
combination showing:

  - Middle axial, coronal, and sagittal slices for the uncorrected image,
    the estimated bias field, and the corrected image (displayed side-by-side
    in columns).
  - Intensity histograms for all three images, with median and 99th-percentile
    markers overlaid.
  - A summary statistics table (min, max, mean, and several percentiles) for
    each of the three images.

The PNG is written to the ``qc`` datatype directory.  It is subsequently
collected by ``qc_n4_tune_report`` into an HTML report.

This is a Snakemake script; the ``snakemake`` object is automatically provided
when executed as part of a Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np


def _load_vol(path):
    """Load a NIfTI file and return a 3-D float32 array (channel dims squeezed)."""
    img = nib.load(path)
    data = np.asarray(img.dataobj).astype(np.float32)
    while data.ndim > 3:
        data = data[..., 0]
    return data


def _mid_slices(vol):
    """Return (axial, coronal, sagittal) middle slices of a 3-D volume."""
    z, y, x = vol.shape
    return vol[z // 2, :, :], vol[:, y // 2, :], vol[:, :, x // 2]


def _norm(arr, lo, hi):
    """Linearly normalise *arr* to [0, 1] using *lo* / *hi* bounds."""
    if hi > lo:
        return np.clip((arr.astype(np.float32) - lo) / (hi - lo), 0.0, 1.0)
    return np.zeros_like(arr, dtype=np.float32)


def _compute_stats(vol):
    """Return a dict of summary statistics for a 3-D volume."""
    flat = vol.ravel()
    pctls = np.percentile(flat, [1, 5, 25, 50, 75, 95, 99])
    return {
        "min": float(flat.min()),
        "max": float(flat.max()),
        "mean": float(flat.mean()),
        "p1": float(pctls[0]),
        "p5": float(pctls[1]),
        "p25": float(pctls[2]),
        "p50": float(pctls[3]),
        "p75": float(pctls[4]),
        "p95": float(pctls[5]),
        "p99": float(pctls[6]),
    }


def main():
    n4_spline = snakemake.params.n4_spline
    n4_iters = snakemake.params.n4_iters
    stain = snakemake.wildcards.stain
    subject = snakemake.wildcards.subject

    print(f"Loading images for N4 tuning QC (spline={n4_spline}, iters={n4_iters})…")
    raw = _load_vol(snakemake.input.nii)
    field = _load_vol(snakemake.input.biasfield)
    corrected = _load_vol(snakemake.input.corrected)

    # Display intensity bounds
    raw_lo = float(np.percentile(raw, 0.5))
    raw_hi = float(np.percentile(raw, 99.5))
    field_lo = float(np.percentile(field, 0.5))
    field_hi = float(np.percentile(field, 99.5))

    # Per-image metadata: (volume, display_lo, display_hi, cmap, label)
    images = [
        (raw, raw_lo, raw_hi, "gray", "Uncorrected"),
        (field, field_lo, field_hi, "RdBu_r", "Bias Field"),
        (corrected, raw_lo, raw_hi, "gray", "Corrected"),
    ]
    orient_labels = ["Axial (mid-Z)", "Coronal (mid-Y)", "Sagittal (mid-X)"]

    stats = {label: _compute_stats(vol) for vol, *_, label in images}

    # -----------------------------------------------------------------------
    # Build figure: 4 rows × 3 columns
    #   rows 0-2 : axial / coronal / sagittal slice montages
    #   row  3   : histogram + in-panel stats text
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(4, 3, figsize=(15, 18))
    fig.suptitle(
        f"N4 Bias Field Correction — Parameter Tuning QC\n"
        f"Subject: {subject}  |  Stain: {stain}  |  "
        f"Spline spacing: {n4_spline} mm  |  Iterations per level: {n4_iters}",
        fontsize=12,
        fontweight="bold",
        y=0.99,
    )

    # Slice panels (rows 0-2)
    for col, (vol, lo, hi, cmap, label) in enumerate(images):
        slices = _mid_slices(vol)
        for row, (sl, orient) in enumerate(zip(slices, orient_labels)):
            ax = axes[row, col]
            ax.imshow(_norm(sl, lo, hi), cmap=cmap, origin="lower", aspect="auto")
            ax.set_title(f"{label}\n{orient}", fontsize=9)
            ax.axis("off")

    # Histogram + stats panels (row 3)
    hist_colors = ["steelblue", "darkorange", "forestgreen"]
    for col, (vol, lo, hi, _, label) in enumerate(images):
        ax = axes[3, col]
        flat = vol.ravel()
        hist_lo = max(float(flat.min()), lo - abs(hi - lo) * 0.05)
        hist_hi = hi + abs(hi - lo) * 0.05
        counts, edges = np.histogram(flat, bins=200, range=(hist_lo, hist_hi))
        centers = (edges[:-1] + edges[1:]) / 2
        ax.bar(
            centers,
            counts,
            width=edges[1] - edges[0],
            color=hist_colors[col],
            alpha=0.75,
        )
        s = stats[label]
        ax.axvline(
            s["p50"],
            color="red",
            linestyle="--",
            lw=1.2,
            label=f"p50={s['p50']:.2f}",
        )
        ax.axvline(
            s["p99"],
            color="purple",
            linestyle="--",
            lw=1.2,
            label=f"p99={s['p99']:.2f}",
        )
        ax.legend(fontsize=7, loc="upper right")
        ax.set_xlabel("Intensity", fontsize=8)
        ax.set_ylabel("Voxel count", fontsize=8)
        ax.set_title(f"{label} — histogram", fontsize=9)
        ax.tick_params(labelsize=7)

        # Stats text below histogram
        stat_lines = (
            f"min={s['min']:.2f}  max={s['max']:.2f}  mean={s['mean']:.2f}\n"
            f"p1={s['p1']:.2f}  p5={s['p5']:.2f}  p25={s['p25']:.2f}\n"
            f"p50={s['p50']:.2f}  p75={s['p75']:.2f}  p95={s['p95']:.2f}  p99={s['p99']:.2f}"
        )
        ax.text(
            0.01,
            -0.32,
            stat_lines,
            transform=ax.transAxes,
            fontsize=7,
            verticalalignment="top",
            fontfamily="monospace",
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8),
        )

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(snakemake.output.png, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Saved N4 tuning QC PNG to {snakemake.output.png}")


if __name__ == "__main__":
    main()
