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

import dask
from dask.diagnostics import ProgressBar
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

# Get wildcards
atlas_seg = snakemake.wildcards.seg

# Get level from wildcards if available
level = int(snakemake.wildcards.get("level", 0))

# Set up dask for parallel processing
dask.config.set(scheduler="threads", num_workers=snakemake.threads)

# Create output directory
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Load the atlas with labels
atlas = ZarrNiiAtlas.from_files(
    input_dseg,
    input_tsv,
    **{k: v for k, v in zarrnii_kwargs.items() if v is not None},
)

# Load the image data
# Check if input is ome.zarr format or nifti
image = ZarrNii.from_ome_zarr(
    input_zarr,
    level=level,
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
                clean_abbrev = re.sub(r"[^a-zA-Z0-9]+", "_", str(label_abbrev)).strip(
                    "_"
                )
                # Fallback to label index if abbreviation would be empty
                if not clean_abbrev:
                    clean_abbrev = f"idx{label_idx}"
                out_file = Path(output_dir) / (
                    f"seg-{atlas_seg}_label-{clean_abbrev}_patch-{i:04d}.nii"
                )
                patch.to_nifti(str(out_file))

        except ValueError as e:
            # ValueError from sample_region_patches when region has no voxels
            logging.warning(
                f"Skipping label {label_abbrev} (index {label_idx}): "
                f"no voxels in region - {e}"
            )
            continue
