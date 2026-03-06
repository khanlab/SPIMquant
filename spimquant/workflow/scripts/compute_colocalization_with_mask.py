"""Compute colocalization of instance segmentation objects with a binary mask.

This script takes aggregated regionprops from instance segmentation stains (in
template space, with subject-space coordinates retained) and determines whether
each object's centroid falls inside a binary segmentation mask from a mask stain
(e.g. a vessel mask).

The colocalization is performed in subject space, since the mask is stored in
subject-space OME-Zarr format. The subject-space coordinates (pos_x, pos_y,
pos_z) from the regionprops are used to query the mask at each object's centroid.

Input:
    - regionprops_parquet: Aggregated regionprops with subject-space coordinates
      (pos_x, pos_y, pos_z) and template-space coordinates (template_x/y/z).
    - mask: OME-Zarr directory containing the binary segmentation mask for the
      mask stain, in subject space.

Parameters:
    - coord_column_names: Column names for subject-space physical coordinates in
      the regionprops table (as defined in config, e.g. ['pos_x', 'pos_y', 'pos_z']).
      These are used to look up mask values and must match the coordinate system of
      the mask's ZarrNii affine.
    - template_coord_column_names: Column names for template-space coordinates
      (e.g. ['template_x', 'template_y', 'template_z']).
    - mask_stain: Name of the mask stain being colocalized against.
    - zarrnii_kwargs: Additional kwargs passed to ZarrNii.from_ome_zarr.

Output:
    - maskcoloc_parquet: Annotated regionprops table with all original columns
      plus:
        * in_mask (int): 1 if the object centroid is inside the mask, 0 otherwise.
        * mask_stain (str): Name of the mask stain used for colocalization.
        * template_coloc_x/y/z (float): Template-space coordinates for downstream
          atlas-based quantification (same as template_x/y/z).

The output can be used for downstream atlas-based quantification to determine
how many instance objects (e.g. amyloid plaques) are located inside mask regions
(e.g. blood vessels).
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interpn
from zarrnii import ZarrNii

# Load the aggregated instance regionprops (template space, with subject coords)
df = pd.read_parquet(snakemake.input.regionprops_parquet)

# Get coordinate column names
coord_cols = snakemake.params.coord_column_names  # subject-space physical coords
template_coord_cols = snakemake.params.template_coord_column_names  # template-space coords
mask_stain = snakemake.params.mask_stain

# Validate coordinate columns exist
for col in coord_cols + template_coord_cols:
    if col not in df.columns:
        raise ValueError(
            f"Coordinate column '{col}' not found in regionprops data. "
            f"Available columns: {df.columns.tolist()}"
        )

# Load the mask image (subject space, binary mask)
mask_znimg = ZarrNii.from_ome_zarr(
    snakemake.input.mask,
    level=0,
    **snakemake.params.zarrnii_kwargs,
)

# Get the inverse affine to convert subject-space physical coordinates to voxel indices.
# The ZarrNii affine maps from voxel indices (ZYX order) to physical space.
# We invert it to go from physical coordinates → voxel indices (ZYX order).
ras_to_vox = mask_znimg.affine.invert()

# Extract subject-space physical coordinates from the regionprops (shape: N x 3).
# The coord_cols are ordered as defined in the config (pos_x, pos_y, pos_z),
# corresponding to the physical coordinate axes of the image.
points_phys = df[coord_cols].values  # N x 3

# Convert physical coordinates → voxel indices (ZYX order).
# apply_transform expects (3, N) input and returns (3, N) output.
# The ZarrNii affine maps ZYX voxels → physical space, so its inverse maps
# physical space → ZYX voxels. The output rows correspond to [Z, Y, X] voxel indices.
points_vox = ras_to_vox.apply_transform(points_phys.T)  # shape: (3, N), ZYX voxel coords

# Determine in_mask values for each point using nearest-neighbor lookup
# (appropriate for binary masks)

# Initialize in_mask array with zeros (default: outside mask)
in_mask = np.zeros(len(df), dtype=np.int8)

if len(df) > 0:
    # Use get_bounded_subregion and interpn for efficient lookup.
    # get_bounded_subregion expects (3, N) or (4, N) input in voxel (ZYX) order.
    # It loads only the subregion of the mask containing the query points.
    grid_points, subvol = mask_znimg.get_bounded_subregion(points_vox)

    if subvol is not None:
        # subvol shape: (C, Z, Y, X) or (Z, Y, X) depending on mask format
        # Squeeze channel dimension if present
        if subvol.ndim == 4:
            subvol = subvol.squeeze(0)

        # Interpolate mask values at each point using nearest-neighbor method.
        # interpn expects:
        #   - grid_points: tuple of 1D arrays (one per dimension, in ZYX order)
        #   - values: N-D array in ZYX order
        #   - xi: query points as (N, 3) array in ZYX order
        mask_values = interpn(
            grid_points,
            subvol.astype(np.float32),
            points_vox[:3, :].T,  # (N, 3) in ZYX order, matching grid_points
            method="nearest",
            bounds_error=False,
            fill_value=0.0,
        )

        # Threshold to binary (handle any float values from the mask)
        in_mask = (mask_values > 0.5).astype(np.int8)

# Annotate the regionprops dataframe with mask colocalization info
df = df.copy()
df["in_mask"] = in_mask
df["mask_stain"] = mask_stain

# Add template_coloc coordinates (same as the instance object's template coordinates)
# These are used by downstream atlas mapping rules
template_coloc_col_names = ["template_coloc_x", "template_coloc_y", "template_coloc_z"]
for coloc_col, tpl_col in zip(template_coloc_col_names, template_coord_cols):
    df[coloc_col] = df[tpl_col]

# Save the annotated regionprops
df.to_parquet(snakemake.output.maskcoloc_parquet, index=False)

n_in_mask = int(in_mask.sum())
print("Mask colocalization analysis complete:")
print(f"  Total instance objects: {len(df)}")
print(f"  Mask stain: {mask_stain}")
print(f"  Objects inside mask: {n_in_mask} ({100 * n_in_mask / max(len(df), 1):.1f}%)")
print(f"  Objects outside mask: {len(df) - n_in_mask}")
if "stain" in df.columns:
    for stain, grp in df.groupby("stain"):
        n_in = int(grp["in_mask"].sum())
        print(
            f"    {stain}: {n_in}/{len(grp)} ({100 * n_in / max(len(grp), 1):.1f}%) inside mask"
        )
