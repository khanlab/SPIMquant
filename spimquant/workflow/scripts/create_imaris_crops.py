"""Create Imaris datasets from zarr data based on atlas region bounding boxes.

This script extracts crops from SPIM zarr data based on bounding boxes of
specified atlas regions. It uses the zarrnii library to load the atlas with
its label lookup table, get bounding boxes for specified regions, crop the
image data, and save as Imaris datasets.

The output is a directory containing Imaris datasets named:
    seg-{atlas_seg}_label-{labelabbrev}.ims
"""

import logging
import re
from pathlib import Path
import numpy as np

import dask
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii, ZarrNiiAtlas

# Set up logging for snakemake scripts
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

# Get input from spim
input_zarr = snakemake.input.spim

input_dseg = snakemake.input.dseg
input_tsv = snakemake.input.label_tsv
output_dir = snakemake.output.crops_dir

# Parameters
crop_labels = snakemake.params.crop_labels
zarrnii_kwargs = snakemake.params.zarrnii_kwargs

# Get wildcards
atlas_seg = snakemake.wildcards.seg

# Get level from wildcards if available
target_level = int(snakemake.wildcards.level)
hires_level = int(snakemake.params.hires_level)

downsampling_level = target_level - hires_level
if downsampling_level < 0:
    raise ValueError(
        "Target level for create_imaris_crops is smaller than the input level!"
    )


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
image = ZarrNii.from_ome_zarr(
    input_zarr,
    level=downsampling_level,
    **{k: v for k, v in zarrnii_kwargs.items() if v is not None},
)

# Determine which labels to use for crops
if crop_labels is None:
    # Use all non-background labels from the atlas
    labels_to_use = atlas.labels_df[atlas.labels_df[atlas.label_column] > 0][
        [atlas.label_column, atlas.abbrev_column]
    ].values.tolist()
else:
    # Use specified labels - can be indices or names/abbreviations
    labels_to_use = []
    for label in crop_labels:
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

# Extract crops for each label
with ProgressBar():
    for label_idx, label_abbrev in labels_to_use:
        try:
            # Get bounding box for this region
            bbox_min, bbox_max = atlas.get_region_bounding_box(
                region_ids=int(label_idx)
            )
            # Crop the image using the bounding box
            cropped = image.crop(bbox_min, bbox_max, physical_coords=True)

            logging.info(f"cropped shape for {label_abbrev} is {cropped.shape}")
            if any(d > 5000 for d in cropped.shape):
                raise ValueError(
                    f"Cropped image too large, shape={cropped.shape}, skipping"
                )

            # Clean label namn for filename: remove non-alphanumeric chars
            clean_abbrev = re.sub(r"[^a-zA-Z0-9]+", "", str(label_abbrev))
            # Fallback to label index if name would be empty
            if not clean_abbrev:
                clean_abbrev = f"idx{label_idx}"
            subject = snakemake.wildcards.subject
            out_file = Path(output_dir) / (
                f"sub-{subject}_seg-{atlas_seg}_label-{clean_abbrev}_SPIM.ims"
            )
            
            # Save as Imaris dataset
            cropped.to_imaris(str(out_file))

            logging.info(
                f"Created Imaris crop for label {label_abbrev} (index {label_idx}): {out_file}"
            )

        except ValueError as e:
            # ValueError from get_region_bounding_box when region has no voxels
            logging.warning(
                f"Skipping label {label_abbrev} (index {label_idx}): "
                f"no voxels in region - {e}"
            )
            continue
        except Exception as e:
            # Catch any other errors
            logging.error(
                f"Error processing label {label_abbrev} (index {label_idx}): {e}"
            )
            continue
