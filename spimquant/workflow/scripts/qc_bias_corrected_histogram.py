"""Bias-corrected per-channel intensity histogram QC for SPIM data.

Samples random full-resolution patches from within the brain mask, applies
bias field correction patch-wise by upsampling the downsampled correction map,
and generates a four-panel intensity histogram of the corrected intensities.

This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import logging

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import zoom

from dask_setup import get_dask_client
from zarrnii import ZarrNii, ZarrNiiAtlas

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# Integer label used in the single-label brain-mask atlas.
BRAIN_REGION_ID = 1

# Multiplicative headroom added to axis limits when computing display bounds.
DISPLAY_MARGIN_FACTOR = 1.05

# Minimum value used as denominator floor when dividing by the bias field,
# preventing division by zero.
BIASFIELD_FLOOR = np.finfo(np.float32).eps


def main():
    stain = snakemake.wildcards.stain
    n_patches = snakemake.params.n_patches
    patch_size = tuple(snakemake.params.patch_size)
    seed = snakemake.params.seed
    hist_bins = snakemake.params.hist_bins
    hist_range = snakemake.params.hist_range
    biasfield_zarr_level = snakemake.params.biasfield_zarr_level
    zarrnii_kwargs = {
        k: v for k, v in snakemake.params.zarrnii_kwargs.items() if v is not None
    }

    # Patch size for the biasfield at its downsampled level.
    # Since each pyramid level halves the voxel count per axis, a patch of
    # `patch_size` voxels at level 0 corresponds to
    # `patch_size / 2**biasfield_zarr_level` voxels at the downsampled level,
    # covering the same physical extent.
    biasfield_patch_size = tuple(
        max(1, p // (2**biasfield_zarr_level)) for p in patch_size
    )

    with get_dask_client("threads", snakemake.threads):
        # Load brain mask as a ZarrNiiAtlas with a single "brain" label.
        brain_znii = ZarrNii.from_file(
            snakemake.input.brain_mask, **zarrnii_kwargs
        )
        labels_df = pd.DataFrame(
            {"index": [BRAIN_REGION_ID], "name": ["brain"], "abbreviation": ["brain"]}
        )
        atlas = ZarrNiiAtlas.create_from_dseg(brain_znii, labels_df)

        # Sample patch centers uniformly within the brain mask (physical coords).
        logging.info(f"Sampling {n_patches} patch centers from brain mask ...")
        centers = atlas.sample_region_patches(
            n_patches=n_patches,
            region_ids=BRAIN_REGION_ID,
            seed=seed,
        )
        logging.info(f"Sampled {len(centers)} centers.")

        # Load raw SPIM at level 0 (full resolution) for patch extraction.
        znimg_raw = ZarrNii.from_ome_zarr(
            snakemake.input.spim,
            level=0,
            channel_labels=[stain],
            **zarrnii_kwargs,
        )

        # Load biasfield at a downsampled pyramid level within the biasfield zarr.
        znimg_biasfield = ZarrNii.from_ome_zarr(
            snakemake.input.biasfield,
            level=biasfield_zarr_level,
        )

        # Collect corrected intensities patch by patch.
        all_intensities = []

        for i, center in enumerate(centers):
            try:
                # Extract full-resolution raw patch.
                # crop_centered may return a single ZarrNii or a list; normalise
                # to list so the rest of the loop is consistent.
                raw_patch = znimg_raw.crop_centered([center], patch_size=patch_size)
                if not isinstance(raw_patch, list):
                    raw_patch = [raw_patch]
                raw_np = np.squeeze(raw_patch[0].data.compute()).astype(np.float32)

                # Extract corresponding biasfield patch at the downsampled level.
                # Using the same physical center but a proportionally smaller voxel
                # count so both patches cover the same physical extent.
                bf_patch = znimg_biasfield.crop_centered(
                    [center], patch_size=biasfield_patch_size
                )
                if not isinstance(bf_patch, list):
                    bf_patch = [bf_patch]
                bf_np = np.squeeze(bf_patch[0].data.compute()).astype(np.float32)

                # Upsample the biasfield patch to match the raw patch spatial shape.
                # Linear interpolation (order=1) is chosen because the bias field
                # is a smooth low-frequency signal; linear upsampling avoids the
                # ringing artefacts that higher-order splines can introduce.
                if raw_np.shape != bf_np.shape:
                    zoom_factors = tuple(
                        r / b for r, b in zip(raw_np.shape, bf_np.shape)
                    )
                    bf_np = zoom(bf_np, zoom_factors, order=1)

                # Apply bias field correction: corrected = raw / biasfield.
                corrected = raw_np / np.maximum(bf_np, BIASFIELD_FLOOR)
                all_intensities.append(corrected.ravel())

                logging.info(f"Processed patch {i + 1}/{len(centers)}")

            except (ValueError, IndexError, RuntimeError) as e:
                logging.warning(f"Skipping patch {i + 1}: {e}")
                continue

    if not all_intensities:
        logging.warning("No valid patches collected; producing empty histogram.")
        all_intensities = [np.zeros(1, dtype=np.float32)]

    intensities = np.concatenate(all_intensities)

    # Compute histogram from the collected corrected intensities.
    hist_counts, bin_edges = np.histogram(
        intensities, bins=hist_bins, range=tuple(hist_range)
    )
    hist_counts = hist_counts.astype(float)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = bin_edges[1] - bin_edges[0]

    total_voxels = hist_counts.sum()
    max_range = hist_range[1]

    nonzero_mask = hist_counts > 0
    disp_max = (
        float(bin_centers[nonzero_mask][-1]) * DISPLAY_MARGIN_FACTOR
        if nonzero_mask.any()
        else max_range
    )
    sat_fraction = (
        float(hist_counts[-1]) / total_voxels * 100 if total_voxels > 0 else 0.0
    )

    if total_voxels > 0:
        mean_val = float(np.sum(bin_centers * hist_counts) / total_voxels)
        cumsum_norm = np.cumsum(hist_counts) / total_voxels
        p50_val = float(
            bin_centers[min(np.searchsorted(cumsum_norm, 0.50), len(bin_centers) - 1)]
        )
        p99_val = float(
            bin_centers[min(np.searchsorted(cumsum_norm, 0.99), len(bin_centers) - 1)]
        )
    else:
        mean_val = p50_val = p99_val = 0.0

    lin_xlim = p99_val * DISPLAY_MARGIN_FACTOR if total_voxels > 0 else max_range
    visible = hist_counts[bin_centers <= lin_xlim]
    lin_ylim = (
        float(visible.max()) * DISPLAY_MARGIN_FACTOR
        if visible.size and visible.max() > 0
        else 1.0
    )

    subject = snakemake.wildcards.subject

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(
        f"Bias-Corrected Intensity Histogram QC\n"
        f"Subject: {subject}  |  Stain: {stain}  |  "
        f"Patches: {len(centers)}  |  Correction: {snakemake.params.correction_method}",
        fontsize=12,
        fontweight="bold",
    )

    # Panel 1: linear-scale histogram
    ax = axes[0, 0]
    ax.bar(bin_centers, hist_counts, width=bin_width, color="steelblue", alpha=0.75)
    ax.set_xlabel("Corrected intensity")
    ax.set_ylabel("Voxel count")
    ax.set_title("Linear-scale histogram")
    ax.set_xlim(0, lin_xlim)
    ax.set_ylim(0, lin_ylim)

    # Panel 2: log-scale histogram
    ax = axes[0, 1]
    log_counts = np.where(hist_counts > 0, np.log10(hist_counts), np.nan)
    ax.bar(bin_centers, log_counts, width=bin_width, color="darkorange", alpha=0.75)
    ax.set_xlabel("Corrected intensity")
    ax.set_ylabel("log\u2081\u2080(voxel count)")
    ax.set_title("Log-scale histogram")
    ax.set_xlim(0, disp_max)

    # Panel 3: cumulative distribution
    ax = axes[1, 0]
    if total_voxels > 0:
        cumsum_pct = cumsum_norm * 100
        ax.plot(bin_centers, cumsum_pct, color="forestgreen", lw=1.5)
        ax.axvline(
            x=p50_val,
            color="purple",
            linestyle="--",
            alpha=0.7,
            label=f"Median ({p50_val:.1f})",
        )
        ax.axvline(
            x=p99_val,
            color="red",
            linestyle="--",
            alpha=0.7,
            label=f"99th pctile ({p99_val:.1f})",
        )
        ax.legend(fontsize=8)
    ax.set_xlabel("Corrected intensity")
    ax.set_ylabel("Cumulative voxels (%)")
    ax.set_title("Cumulative distribution")
    ax.set_ylim(0, 105)
    ax.set_xlim(0, disp_max)

    # Panel 4: summary statistics
    ax = axes[1, 1]
    ax.axis("off")
    summary_text = (
        f"Sampled voxels:    {int(total_voxels):>14,}\n"
        f"Patches:           {len(centers):>14,}\n"
        f"Mean intensity:    {mean_val:>14.2f}\n"
        f"Median (50th):     {p50_val:>14.2f}\n"
        f"99th percentile:   {p99_val:>14.2f}\n"
        f"Max range:         {max_range:>14.1f}\n"
        f"Saturation frac.:  {sat_fraction:>13.3f}%"
    )
    ax.text(
        0.1,
        0.55,
        summary_text,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="center",
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8),
    )
    ax.set_title("Summary statistics")

    plt.tight_layout()
    plt.savefig(snakemake.output.png, dpi=150, bbox_inches="tight")
    plt.close()
    logging.info(f"Saved bias-corrected histogram QC to {snakemake.output.png}")


if __name__ == "__main__":
    main()
