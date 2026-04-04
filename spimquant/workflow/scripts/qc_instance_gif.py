"""Instance-centred animated GIF QC.

For every detected instance (plaque, cell, colocalized pair, …) crops the
raw SPIM image around that instance and assembles the crops into two
animated GIFs:

* **random** – instances in a randomised order so rare large/small objects
  are not systematically relegated to the end of the animation.
* **sorted** – instances sorted in ascending order of equivalent spherical
  radius so the viewer can sweep from small to large objects.

Each GIF frame shows the crop for **all** stain channels placed side-by-side
(one subplot per channel).  A hollow circle whose diameter matches the
equivalent spherical diameter of the instance is drawn at the centre of
every subplot.  A text annotation below the subplots records the instance
index, its subject-space position, the atlas ROI it falls in, and its
equivalent radius.

Inputs are controlled by ``snakemake.params.instance_type``:

* ``"stain"`` – per-stain instances read from the aggregated regionprops
  parquet (columns ``pos_x``, ``pos_y``, ``pos_z``, ``nvoxels``).  Atlas
  labels are looked up at runtime from the atlas dseg and label TSV.
* ``"coloc"`` – colocalized pairs read from the coloc parquet (columns
  ``pos_coloc_x``, ``pos_coloc_y``, ``pos_coloc_z``, ``radius_a``,
  ``radius_b``).  Atlas labels are likewise looked up at runtime.

This is a Snakemake script that expects the ``snakemake`` object to be
available.
"""

from __future__ import annotations

import io
import random

import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from zarrnii import ZarrNii, ZarrNiiAtlas

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _estimate_global_percentiles(arr, pct_low=1, pct_high=99):
    """Return (lo, hi) percentile values for intensity normalisation."""
    lo = float(np.percentile(arr[arr > 0], pct_low)) if np.any(arr > 0) else 0.0
    hi = float(np.percentile(arr, pct_high))
    return lo, hi


def _apply_norm(arr, lo, hi):
    """Clip and rescale *arr* to [0, 1]."""
    if hi > lo:
        return np.clip((arr.astype(np.float32) - lo) / (hi - lo), 0.0, 1.0)
    return np.zeros_like(arr, dtype=np.float32)


def _load_channel_images(spim_path, channels, level):
    """Return a dict {channel: ZarrNii} for all requested channels."""
    imgs = {}
    for ch in channels:
        imgs[ch] = ZarrNii.from_ome_zarr(
            spim_path,
            level=level,
            channel_labels=[ch],
        )
    return imgs


def _global_norms(imgs, channels):
    """Compute per-channel (lo, hi) from a coarse level for normalisation."""
    norms = {}
    for ch in channels:
        # load a very coarse version for statistics
        try:
            coarse = ZarrNii.from_ome_zarr(
                snakemake.input.spim,
                level=(int(snakemake.params.level) + 5),
                channel_labels=[ch],
            )
            arr = coarse.data.compute()
        except Exception:
            arr = imgs[ch].data.compute()
        lo, hi = _estimate_global_percentiles(arr)
        norms[ch] = (lo, hi)
    return norms


def _add_atlas_labels(df, pos_cols, dseg_nii, label_tsv):
    """Use ZarrNiiAtlas to add 'atlas_index' and 'atlas_name' columns."""
    atlas = ZarrNiiAtlas.from_files(dseg_nii, label_tsv)
    regionprops_dict = df[pos_cols].copy()
    regionprops_dict = regionprops_dict.rename(
        columns={pos_cols[0]: "pos_x", pos_cols[1]: "pos_y", pos_cols[2]: "pos_z"}
    ).to_dict(orient="list")

    # label_region_properties needs at least one numeric column alongside coords
    # We add a dummy column to satisfy the interface if needed
    regionprops_dict["_dummy"] = [0] * len(df)

    try:
        df_labeled, _ = atlas.label_region_properties(
            regionprops_dict,
            coord_column_names=["pos_x", "pos_y", "pos_z"],
            include_names=True,
        )
        # Find the atlas index/name columns added by label_region_properties
        for candidate in ("index", "label_index", "atlas_index", "region_index"):
            if candidate in df_labeled.columns:
                df["atlas_index"] = df_labeled[candidate].values
                break
        for candidate in ("name", "label_name", "atlas_name", "region_name"):
            if candidate in df_labeled.columns:
                df["atlas_name"] = df_labeled[candidate].values
                break
    except Exception:
        pass  # atlas labeling is best-effort; we just won't have names

    if "atlas_index" not in df.columns:
        df["atlas_index"] = np.nan
    if "atlas_name" not in df.columns:
        df["atlas_name"] = "?"
    return df


def _draw_hollow_circle(ax, cx, cy, radius_px, color="yellow", linewidth=1.5):
    """Draw a hollow circle on *ax* at (cx, cy) with the given pixel radius."""
    if radius_px > 0:
        circle = plt.Circle(
            (cx, cy),
            radius=radius_px,
            fill=False,
            edgecolor=color,
            linewidth=linewidth,
        )
        ax.add_patch(circle)


def _make_frame(row, imgs, channels, norms, pos_cols, patch_size, instance_type):
    """Render one GIF frame as a PIL Image."""
    pos = tuple(float(row[c]) for c in pos_cols)

    # Equivalent radius in voxels
    if instance_type == "coloc":
        radius_vox = float((row.get("radius_a", 1) + row.get("radius_b", 1)) / 2.0)
    else:
        nvox = float(row.get("nvoxels", 1))
        radius_vox = (3.0 * nvox / (4.0 * np.pi)) ** (1.0 / 3.0)

    n_ch = len(channels)
    fig, axes = plt.subplots(
        1, n_ch, figsize=(n_ch * 2.8, 3.5), constrained_layout=False
    )
    if n_ch == 1:
        axes = [axes]

    for ax, ch in zip(axes, channels):
        try:
            crop = imgs[ch].crop_centered(pos, patch_size=(patch_size, patch_size, 1))
            sl = crop.data[0, :, :].squeeze().compute()
        except Exception:
            sl = np.zeros((patch_size, patch_size), dtype=np.float32)

        lo, hi = norms.get(ch, (0.0, 1.0))
        sl_norm = _apply_norm(sl, lo, hi)

        ax.imshow(sl_norm, cmap="gray", interpolation="nearest")
        ax.set_title(ch, fontsize=7, pad=2)
        ax.set_xticks([])
        ax.set_yticks([])

        H, W = sl_norm.shape
        cy, cx = H / 2.0, W / 2.0
        _draw_hollow_circle(ax, cx, cy, radius_vox)

    # Annotation text below the subplots
    atlas_name = str(row.get("atlas_name", "?"))
    atlas_idx = row.get("atlas_index", "")
    roi_str = (
        f"{atlas_name} ({atlas_idx})"
        if atlas_idx not in ("", "?", np.nan)
        else atlas_name
    )
    p0, p1, p2 = pos
    ann = (
        f"idx={int(row.get('_row_idx', 0))}  |  "
        f"pos=({p0:.1f},{p1:.1f},{p2:.1f})  |  "
        f"roi={roi_str}  |  "
        f"r={radius_vox:.1f}vox"
    )
    fig.text(
        0.5,
        0.01,
        ann,
        ha="center",
        va="bottom",
        fontsize=6,
        family="monospace",
        wrap=True,
    )

    fig.subplots_adjust(bottom=0.12, top=0.92, left=0.02, right=0.98, wspace=0.05)

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=72)
    plt.close(fig)
    buf.seek(0)
    pil_img = Image.open(buf).copy()
    buf.close()
    return pil_img


def _save_gif(frames, output_path, duration_ms=600):
    """Save a list of PIL Images as an animated GIF."""
    if not frames:
        # Write a blank placeholder
        blank = Image.new("RGB", (200, 100), color="white")
        blank.save(output_path, format="GIF")
        return

    # Convert to palette mode for smaller file size
    palette_frames = []
    for f in frames:
        palette_frames.append(f.convert("P", palette=Image.ADAPTIVE, colors=256))

    palette_frames[0].save(
        output_path,
        format="GIF",
        save_all=True,
        append_images=palette_frames[1:],
        loop=0,
        duration=duration_ms,
        optimize=False,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    instance_type = snakemake.params.instance_type  # "stain" or "coloc"
    channels = list(snakemake.params.channels)
    patch_size = int(snakemake.params.patch_size)
    max_instances = int(snakemake.params.max_instances)
    level = int(snakemake.params.level)
    seed = int(snakemake.params.get("seed", 42))

    stain_wc = getattr(snakemake.wildcards, "stain", instance_type)
    desc_wc = snakemake.wildcards.desc
    subject_wc = snakemake.wildcards.subject

    # ------------------------------------------------------------------
    # Load instance table
    # ------------------------------------------------------------------
    if instance_type == "coloc":
        df = pd.read_parquet(snakemake.input.instance_parquet)
        pos_cols = ["pos_coloc_x", "pos_coloc_y", "pos_coloc_z"]
    else:
        df = pd.read_parquet(snakemake.input.instance_parquet)
        # Filter to the requested stain
        if "stain" in df.columns:
            df = df[df["stain"] == stain_wc].copy()
        pos_cols = ["pos_x", "pos_y", "pos_z"]

    # Drop rows with missing positions
    df = df.dropna(subset=pos_cols).reset_index(drop=True)

    # Add atlas labels
    df = _add_atlas_labels(
        df, pos_cols, snakemake.input.dseg_nii, snakemake.input.label_tsv
    )

    n_instances = len(df)

    if n_instances == 0:
        # Write placeholder GIFs
        blank = Image.new("RGB", (300, 100), color="white")
        for out in (snakemake.output.gif_random, snakemake.output.gif_sorted):
            blank.save(out, format="GIF")
        print(
            f"No instances found for {instance_type} "
            f"(stain={stain_wc}, desc={desc_wc}). "
            "Wrote blank GIFs."
        )
        return

    # Sub-sample if needed (deterministic)
    if n_instances > max_instances:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n_instances, size=max_instances, replace=False)
        df = df.iloc[sorted(idx)].reset_index(drop=True)

    # Preserve original row index for annotation
    df["_row_idx"] = df.index

    # ------------------------------------------------------------------
    # Compute equivalent radius for sorting
    # ------------------------------------------------------------------
    if instance_type == "coloc":
        df["_radius_sort"] = (df.get("radius_a", 1) + df.get("radius_b", 1)) / 2.0
    else:
        nvox = df.get("nvoxels", pd.Series(np.ones(len(df))))
        df["_radius_sort"] = (3.0 * nvox / (4.0 * np.pi)) ** (1.0 / 3.0)

    # ------------------------------------------------------------------
    # Load SPIM channels and compute global normalisations
    # ------------------------------------------------------------------
    imgs = _load_channel_images(snakemake.input.spim, channels, level)
    norms = _global_norms(imgs, channels)

    # ------------------------------------------------------------------
    # Build two orderings
    # ------------------------------------------------------------------
    rng_state = random.Random(seed)
    random_order = list(df.index)
    rng_state.shuffle(random_order)
    df_random = df.loc[random_order].reset_index(drop=True)

    df_sorted = df.sort_values("_radius_sort").reset_index(drop=True)

    # ------------------------------------------------------------------
    # Render frames and save GIFs
    # ------------------------------------------------------------------
    for df_ordered, out_path in (
        (df_random, snakemake.output.gif_random),
        (df_sorted, snakemake.output.gif_sorted),
    ):
        frames = []
        for _, row in df_ordered.iterrows():
            frame = _make_frame(
                row, imgs, channels, norms, pos_cols, patch_size, instance_type
            )
            frames.append(frame)
        _save_gif(frames, out_path)
        print(
            f"Saved {len(frames)}-frame GIF ({instance_type}, {stain_wc}, {desc_wc}, "
            f"subject={subject_wc}) → {out_path}"
        )


if __name__ == "__main__":
    main()
