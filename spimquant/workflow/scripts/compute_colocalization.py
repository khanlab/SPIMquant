"""Compute colocalization between objects from different channels.

This script takes aggregated regionprops from multiple stains (in template space)
and performs KDTree-based spatial colocalization analysis across all channel pairs.
For each pair of nearby objects from different channels, it records their relationship
and computes a colocalization coordinate.

Input:
    - regionprops_aggregated_parquet: Aggregated regionprops with 'stain' column and 
      template space coordinates (template_x, template_y, template_z)

Parameters:
    - search_radius_multiplier: Multiplier for object radius to define search distance
      (default: 3.0). A value of 3.0 means searching within 3x the object's radius.
    - overlap_threshold: Minimum overlap ratio to record a colocalization
      (default: 0.0). Value of 0.0 records all potential overlaps.

Output:
    - coloc_parquet: Table of colocalized object pairs with:
      * object_id_a, object_id_b: IDs of the paired objects
      * stain_a, stain_b: Channels/stains of each object
      * radius_a, radius_b: Estimated radii (from nvoxels)
      * nvoxels_a, nvoxels_b: Volume in voxels
      * distance: Euclidean distance between centroids
      * overlap_ratio: Estimated overlap (1 - distance/sum_radii)
      * template_coloc_x/y/z: Colocalization coordinate (midpoint)

The output can be used for downstream spatial statistics including KDE, 
histograms, voxelization, or further colocalization analysis.
"""

import numpy as np
import pandas as pd
from scipy.spatial import KDTree


# Get configuration parameters with defaults
# Search radius multiplier: determines how far to look for potential colocalizations
# A value of 3.0 is conservative - it searches within 3x the object's radius
# Smaller values (e.g., 1.5) give tighter colocalization, larger values (e.g., 5.0) are more permissive
SEARCH_RADIUS_MULTIPLIER = 1.0

# Overlap threshold: minimum overlap ratio to record a colocalization
# 0.0 records all overlaps (distance < sum of radii)
# Higher values (e.g., 0.5) require more significant overlap
OVERLAP_THRESHOLD = 0.0

# Load the aggregated regionprops data (in template space)
df = pd.read_parquet(snakemake.input.regionprops_parquet)


# Get coordinate column names from config
coord_cols = snakemake.params.coord_column_names
template_coord_cols = snakemake.params.template_coord_column_names

# Validate that we have exactly 3 coordinate columns for 3D analysis
if len(coord_cols) != 3:
    raise ValueError(
        f"Expected exactly 3 coordinate columns for 3D analysis, got {len(coord_cols)}: {coord_cols}"
    )

# Validate that coordinate columns exist
for col in coord_cols:
    if col not in df.columns:
        raise ValueError(
            f"Coordinate column '{col}' not found in regionprops data. Available columns: {df.columns.tolist()}"
        )

# Validate that we have exactly 3 coordinate columns for 3D analysis
if len(template_coord_cols) != 3:
    raise ValueError(
        f"Expected exactly 3 coordinate columns for 3D analysis, got {len(template_coord_cols)}: {template_coord_cols}"
    )

# Validate that coordinate columns exist
for col in template_coord_cols:
    if col not in df.columns:
        raise ValueError(
            f"Coordinate column '{col}' not found in regionprops data. Available columns: {df.columns.tolist()}"
        )


# Validate that required columns exist
required_cols = ["stain", "nvoxels"]
for col in required_cols:
    if col not in df.columns:
        raise ValueError(
            f"Required column '{col}' not found in regionprops data. Available columns: {df.columns.tolist()}"
        )


def estimate_radius_from_nvoxels(nvoxels, voxel_size=1.0):
    """Estimate radius from number of voxels assuming spherical shape.

    Parameters
    ----------
    nvoxels : array-like
        Number of voxels per object
    voxel_size : float, optional
        Size of each voxel in physical units (default: 1.0)

    Returns
    -------
    radius : array-like
        Estimated radius for each object
    """
    # Volume = (4/3) * pi * r^3 = nvoxels * voxel_size^3
    # Solving for r: r = ((3 * nvoxels * voxel_size^3) / (4 * pi))^(1/3)
    volume = nvoxels * (voxel_size**3)
    radius = np.cbrt((3 * volume) / (4 * np.pi))
    return radius


# Add unique object IDs and compute radii
# Use compound IDs to ensure uniqueness across stains
# Use vectorized operation for better performance
df["object_id"] = df["stain"] + "_" + df.index.astype(str)
df["radius"] = estimate_radius_from_nvoxels(df["nvoxels"].values, voxel_size=0.002)


# Get list of unique stains/channels
stains = df["stain"].unique()

# Prepare list to store colocalization results
coloc_results = []


# Perform pairwise colocalization across all channel pairs
for i, stain_a in enumerate(stains):
    for j, stain_b in enumerate(stains):
        # Only process each pair once (and skip self-pairing)
        if i >= j:
            continue

        # Get objects for each stain
        df_a = df[df["stain"] == stain_a].copy()
        df_b = df[df["stain"] == stain_b].copy()

        # Extract coordinates
        coords_a = df_a[coord_cols].values
        coords_b = df_b[coord_cols].values

        # Get template coords (use these when computing final coords)
        template_coords_a = df_a[template_coord_cols].values
        template_coords_b = df_b[template_coord_cols].values

        # Build KDTree for stain_b
        tree_b = KDTree(coords_b)

        # For each object in stain_a, find nearby objects in stain_b
        for idx_a in range(len(df_a)):
            obj_a = df_a.iloc[idx_a]
            pos_a = coords_a[idx_a]
            template_pos_a = template_coords_a[idx_a]

            radius_a = obj_a["radius"]

            # Query the tree for objects within search distance
            # Use configured multiplier to determine search radius
            search_distance = radius_a * SEARCH_RADIUS_MULTIPLIER

            indices = tree_b.query_ball_point(pos_a, r=search_distance)

            # Process each nearby object from stain_b
            for idx_b in indices:
                obj_b = df_b.iloc[idx_b]
                pos_b = coords_b[idx_b]
                template_pos_b = template_coords_b[idx_b]

                radius_b = obj_b["radius"]

                # Calculate distance between centroids
                distance = np.linalg.norm(pos_a - pos_b)

                # Estimate overlap based on distance and radii
                # Overlap is significant if distance < sum of radii
                sum_radii = radius_a + radius_b

                # Handle edge case where both radii are zero
                if sum_radii == 0:
                    overlap_ratio = 0.0
                else:
                    overlap_ratio = max(0, 1 - (distance / sum_radii))

                # Only record if overlap exceeds threshold
                if overlap_ratio > OVERLAP_THRESHOLD:
                    # Calculate colocalization coordinate (midpoint)
                    coloc_coord = (template_pos_a + template_pos_b) / 2.0

                    # Store the colocalization result
                    coloc_results.append(
                        {
                            "object_id_a": obj_a["object_id"],
                            "object_id_b": obj_b["object_id"],
                            "stain_a": stain_a,
                            "stain_b": stain_b,
                            "radius_a": radius_a,
                            "radius_b": radius_b,
                            "nvoxels_a": obj_a["nvoxels"],
                            "nvoxels_b": obj_b["nvoxels"],
                            "distance": distance,
                            "overlapratio": overlap_ratio,
                            "template_coloc_x": coloc_coord[0],
                            "template_coloc_y": coloc_coord[1],
                            "template_coloc_z": coloc_coord[2],
                        }
                    )


# Create output dataframe
if len(coloc_results) > 0:
    df_coloc = pd.DataFrame(coloc_results)
else:
    # Create empty dataframe with correct columns if no colocalizations found
    df_coloc = pd.DataFrame(
        columns=[
            "object_id_a",
            "object_id_b",
            "stain_a",
            "stain_b",
            "radius_a",
            "radius_b",
            "nvoxels_a",
            "nvoxels_b",
            "distance",
            "overlapratio",
            "template_coloc_x",
            "template_coloc_y",
            "template_coloc_z",
        ]
    )

# Save the colocalization results
df_coloc.to_parquet(snakemake.output.coloc_parquet, index=False)

print(f"Colocalization analysis complete:")
print(f"  Total objects: {len(df)}")
print(f"  Stains analyzed: {', '.join(stains)}")
print(f"  Colocalized pairs found: {len(df_coloc)}")
