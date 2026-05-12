"""Create patches from zarr data based on atlas regions.

This script extracts fixed-size patches from zarr data (SPIM or mask) at
locations sampled from specified atlas regions. It uses the zarrnii library
to load the atlas with its label lookup table, sample patch centers from
specified regions, and extract patches from the image data.

The output is a directory containing NIfTI files named:
    seg-{atlas_seg}_label-{labelabbrev}_patch-{patchnum}.nii
"""

import logging
import re
from pathlib import Path

import nibabel as nib
import numpy as np
from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from zarrnii import ZarrNii, ZarrNiiAtlas

# Set up logging for snakemake scripts
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

channel_args = {}
# Determine which input we have (spim, mask, or corrected)
if hasattr(snakemake.input, "spim"):
    input_zarr = snakemake.input.spim
    channel_args = {"channel_labels": [snakemake.wildcards.stain]}
elif hasattr(snakemake.input, "mask"):
    input_zarr = snakemake.input.mask
elif hasattr(snakemake.input, "corrected"):
    input_zarr = snakemake.input.corrected
else:
    raise ValueError("No valid input found (expected spim, mask, or corrected)")

input_dseg = snakemake.input.dseg
input_tsv = snakemake.input.label_tsv
output_dir = snakemake.output.patches_dir

# Parameters
patch_size = tuple(snakemake.params.patch_size)
n_patches = snakemake.params.n_patches
patch_labels = snakemake.params.patch_labels
seed = snakemake.params.seed
zarrnii_kwargs = snakemake.params.zarrnii_kwargs
patch_uint8 = snakemake.params.get("patch_uint8", True)

# Get wildcards
atlas_seg = snakemake.wildcards.seg

# Get level from wildcards if available
target_level = int(snakemake.wildcards.level)
hires_level = int(snakemake.params.hires_level)

downsampling_level = target_level - hires_level
if downsampling_level < 0:
    raise ValueError("Target level for create_patches is smaller than the input level!")

# Create output directory
Path(output_dir).mkdir(parents=True, exist_ok=True)


def convert_nii_to_uint8(
    nii_path, low_pct=0.1, high_pct=99.9, sample_voxels=5_000_000, rng_seed=42
):
    """Convert a NIfTI file to uint8 using robust percentile scaling.

    Scales voxel intensities to the [0, 255] range using the specified
    percentiles for clipping, then overwrites the file in place.
    This avoids endianness issues in downstream applications such as SyGlass.

    Parameters
    ----------
    nii_path : str
        Path to the NIfTI file to convert (overwritten in place).
    low_pct : float
        Lower percentile (0-100) used for intensity clipping. Default is 0.1.
    high_pct : float
        Upper percentile (0-100) used for intensity clipping. Default is 99.9.
    sample_voxels : int
        Maximum number of nonzero voxels to sample for percentile estimation.
        Default is 5,000,000.
    rng_seed : int
        Random seed for reproducible voxel sampling. Default is 42.
    """
    img = nib.load(nii_path)

    # Skip conversion if already uint8
    if img.get_data_dtype() == np.uint8:
        return

    data = np.asanyarray(img.dataobj)

    flat = data.ravel()
    nonzero = flat[flat > 0]

    if nonzero.size == 0:
        out = np.zeros(data.shape, dtype=np.uint8)
        new_header = img.header.copy()
        new_header.set_data_dtype(np.uint8)
        new_header["scl_slope"] = 1
        new_header["scl_inter"] = 0
        nib.save(nib.Nifti1Image(out, img.affine, new_header), nii_path)
        return

    if nonzero.size > sample_voxels:
        rng = np.random.default_rng(rng_seed)
        sample = rng.choice(nonzero, size=sample_voxels, replace=False)
    else:
        sample = nonzero

    lo, hi = np.percentile(sample, [low_pct, high_pct])

    if hi <= lo:
        raise ValueError(f"Bad percentile range: lo={lo}, hi={hi}")

    logging.info(
        "Converting %s to uint8: dtype=%s, %g%%=%g, %g%%=%g",
        nii_path,
        data.dtype,
        low_pct,
        lo,
        high_pct,
        hi,
    )

    data_f = data.astype(np.float32)
    scaled = (data_f - lo) / (hi - lo)
    scaled = np.clip(scaled, 0, 1)
    out = np.round(scaled * 255).astype(np.uint8)

    new_header = img.header.copy()
    new_header.set_data_dtype(np.uint8)
    new_header["scl_slope"] = 1
    new_header["scl_inter"] = 0
    nib.save(nib.Nifti1Image(out, img.affine, new_header), nii_path)


with get_dask_client("threads", snakemake.threads):
    # Load the atlas with labels
    atlas = ZarrNiiAtlas.from_files(
        input_dseg,
        input_tsv,
        **{k: v for k, v in zarrnii_kwargs.items() if v is not None},
    )

    # Load the image data
    # Check if input is ome.zarr format or nifti
    image = ZarrNii.from_file(
        input_zarr,
        level=downsampling_level,
        **channel_args,
        **{k: v for k, v in zarrnii_kwargs.items() if v is not None},
    )

    # Determine which labels to use for patches
    if patch_labels is None:
        # Use all non-background labels from the atlas
        labels_to_use = atlas.labels_df[atlas.labels_df[atlas.label_column] > 0][
            [atlas.label_column, atlas.abbrev_column]
        ].values.tolist()
    else:
        # Use specified labels - can be indices or names/abbreviations
        labels_to_use = []
        for label in patch_labels:
            if isinstance(label, int):
                # Label index provided
                row = atlas.labels_df[atlas.labels_df[atlas.label_column] == label]
                if not row.empty:
                    labels_to_use.append([label, row[atlas.abbrev_column].values[0]])
            else:
                # Label name or abbreviation provided
                row = atlas.labels_df[
                    (atlas.labels_df[atlas.name_column] == label)
                    | (atlas.labels_df[atlas.abbrev_column] == label)
                ]
                if not row.empty:
                    labels_to_use.append(
                        [
                            row[atlas.label_column].values[0],
                            row[atlas.abbrev_column].values[0],
                        ]
                    )

    # Extract patches for each label
    with ProgressBar():
        for label_idx, label_abbrev in labels_to_use:
            try:
                # Sample patch centers from this region
                centers = atlas.sample_region_patches(
                    n_patches=n_patches,
                    region_ids=int(label_idx),
                    seed=seed,
                )

                # Extract patches at each center
                patches = image.crop_centered(centers, patch_size=patch_size)

                # Handle single patch case (when n_patches=1)
                if not isinstance(patches, list):
                    patches = [patches]

                # Save each patch
                for i, patch in enumerate(patches):
                    # Clean label abbreviation for filename: replace non-alphanumeric
                    # chars with underscore, collapse multiple underscores, strip edges
                    clean_abbrev = re.sub(
                        r"[^a-zA-Z0-9]+", "_", str(label_abbrev)
                    ).strip("_")
                    # Fallback to label index if abbreviation would be empty
                    if not clean_abbrev:
                        clean_abbrev = f"idx{label_idx}"
                    out_file = Path(output_dir) / (
                        f"seg-{atlas_seg}_label-{clean_abbrev}_patch-{i:04d}.nii"
                    )
                    patch.to_nifti(str(out_file))
                    if patch_uint8:
                        convert_nii_to_uint8(str(out_file), rng_seed=seed)

            except ValueError as e:
                # ValueError from sample_region_patches when region has no voxels
                logging.warning(
                    f"Skipping label {label_abbrev} (index {label_idx}): "
                    f"no voxels in region - {e}"
                )
                continue
