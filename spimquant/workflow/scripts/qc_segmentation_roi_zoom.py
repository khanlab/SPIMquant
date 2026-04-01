"""ROI-cropped segmentation montage QC.

For each brain region in the atlas parcellation (resampled to subject space),
crops the SPIM image and the segmentation mask to a fixed 2D bounding box 
at the region centroid.   This provides a detail-level view of segmentation 
quality within individual brain regions, complementing the whole-brain 
overview in ``qc_segmentation_overview``.


This is a Snakemake script that expects the ``snakemake`` object to be
available, which is automatically provided when executed as part of a
Snakemake workflow.
"""

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
from zarrnii import ZarrNii, ZarrNiiAtlas
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import zoom


def _percentile_norm(arr, pct_low=1, pct_high=99):
    """Percentile-normalise *arr* to the range [0, 1] using global statistics."""
    lo = np.percentile(arr, pct_low)
    hi = np.percentile(arr, pct_high)
    if hi > lo:
        return np.clip((arr.astype(float) - lo) / (hi - lo), 0.0, 1.0)
    return np.zeros_like(arr, dtype=float)


def _match_shape(source, target_shape, order=1):
    """Zoom *source* array to *target_shape* if shapes differ."""
    if source.shape == target_shape:
        return source
    factors = [t / s for t, s in zip(target_shape, source.shape)]
    return zoom(source, factors, order=order)


def _select_best_z_slice(ff_crop):
    """Return the Z-index with the most field-fraction signal.

    Falls back to the central slice when no signal is present.
    """
    ff_per_z = ff_crop.sum(axis=(0, 1))
    if ff_per_z.max() > 0:
        return int(ff_per_z.argmax())
    return ff_crop.shape[2] // 2


def _get_bounding_box(mask, pad=5):
    """Return index slices for the bounding box of a boolean *mask* with padding.

    Returns ``None`` when the mask is empty.
    """
    indices = np.where(mask)
    if not indices[0].size:
        return None
    shape = mask.shape
    return tuple(
        slice(max(0, int(idx.min()) - pad), min(sz, int(idx.max()) + pad + 1))
        for idx, sz in zip(indices, shape)
    )


def main():
    stain = snakemake.wildcards.stain
    desc = snakemake.wildcards.desc
    subject = snakemake.wildcards.subject
    max_rois = snakemake.params.max_rois
    n_cols = snakemake.params.n_cols

    spim_img = ZarrNii.from_ome_zarr(snakemake.input.spim,level=snakemake.params.level, downsample_near_isotropic=True,channel_labels=[snakemake.wildcards.stain])
    mask_img = ZarrNii.from_ome_zarr(snakemake.input.mask,level=0)
    
    atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg_nii,snakemake.input.label_tsv)

    dseg_data = atlas.dseg.data.compute()


    # Voxel dimensions (mm) for physical aspect-ratio correction - not implemented yet
    # but should be easy with ZarrNii image .scale 
    aspect_axial = 1


    # Global normalisations  - not implemented yet

    # Load atlas label table
    label_df = atlas.labels_df

    # Keep non-background labels that are present in this subject's dseg
    present_ids = set(np.unique(dseg_data)) - {0}
    roi_rows = [
        row
        for _, row in label_df[label_df["index"] > 0].iterrows()
        if int(row["index"]) in present_ids
    ][:max_rois]

    n_rois = len(roi_rows)

    if n_rois == 0:
        fig, ax = plt.subplots(figsize=(18, 12))
        ax.text(
            0.5,
            0.5,
            "No atlas ROIs found in subject",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="gray",
        )
        ax.axis("off")
        plt.savefig(snakemake.output.png, dpi=120, bbox_inches="tight")
        plt.close()
        return

    n_rows = int(np.ceil(n_rois / n_cols))
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3), constrained_layout=True
    )
    fig.suptitle(
        f"ROI Zoom Montage QC\n"
        f"Subject: {subject}  |  Stain: {stain}  |  Method: {desc}",
        fontsize=11,
        fontweight="bold",
    )

    # Normalise axes array to always be 2-D
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes[np.newaxis, :]
    elif n_cols == 1:
        axes = axes[:, np.newaxis]

    for i, row in enumerate(roi_rows):
        ax_row = i // n_cols
        ax_col = i % n_cols
        ax = axes[ax_row, ax_col]

        label_id = int(row["index"])
        label_name = str(row.get("name", label_id))


        #get cropped images for this label
        bbox_min, bbox_max = atlas.get_region_bounding_box(region_ids=label_id)
        center_coord = tuple((x + y) / 2 for x, y in zip(bbox_min, bbox_max))
        spim_crop = spim_img.crop_centered(center_coord, patch_size=(2000,2000,1))
        mask_crop = mask_img.crop_centered(center_coord, patch_size=(2000,2000,1))


        spim_sl = np.rot90(spim_crop.data[0, :, :].squeeze().compute())
        mask_sl = np.rot90(mask_crop.data[0, :, :].squeeze().compute())

        ax.imshow(spim_sl, cmap="gray")
        mask_masked = np.ma.masked_where(mask_sl < 100, mask_sl)
        ax.imshow(mask_masked, cmap="spring", alpha=0.6, vmin=0, vmax=100, aspect=aspect_axial)
        ax.set_title(label_name, fontsize=7, pad=2)
        ax.set_xticks([])
        ax.set_yticks([])

    # Hide unused axes
    for i in range(n_rois, n_rows * n_cols):
        axes[i // n_cols, i % n_cols].axis("off")

    plt.savefig(snakemake.output.png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved ROI zoom montage to {snakemake.output.png}")


if __name__ == "__main__":
    main()
