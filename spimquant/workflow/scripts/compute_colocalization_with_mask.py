"""Compute colocalization of instance-based objects with a mask-based stain.

For each instance in the regionprops (centroid position and nvoxels), this script:

1. Samples the signed distance transform (SDT) of the mask-based stain at the
   instance centroid coordinate using linear interpolation.
2. Calculates an equivalent sphere radius from the instance's nvoxels and the
   physical voxel volume of the SDT image.
3. Adds a new column named ``sdt_{mask_stain}`` = sdt_value - radius:
   - ~0: the instance sphere boundary approximately touches the mask boundary
   - Negative: the instance overlaps with the mask (centroid is inside, or the
     sphere extends into the mask)
   - Positive: the instance is completely outside and away from the mask

The SDT convention used here:
- Negative values: interior of the mask (inside the segmented structure)
- Positive values: exterior of the mask (outside the segmented structure)

Input (via snakemake.input):
    regionprops_parquet: Parquet file with per-instance regionprops including
        centroid columns (pos_x, pos_y, pos_z) and nvoxels.
    dist: OME-Zarr directory containing the signed distance transform of the
        mask-based stain, produced by the ``signed_distance_transform`` rule.

Parameters (via snakemake.params):
    coord_column_names: List of three column names for the x, y, z centroid
        coordinates in the regionprops table (e.g. ['pos_x', 'pos_y', 'pos_z']).
    mask_stain: Name of the mask-based stain, used to name the new column
        (e.g. 'CD31' produces column 'sdt_CD31').
    zarrnii_kwargs: Additional keyword arguments passed to ZarrNii.from_ome_zarr
        (e.g. {'orientation': 'RPI'}).

Output (via snakemake.output):
    maskcoloc_parquet: Parquet file identical to the input regionprops but with
        the additional ``sdt_{mask_stain}`` colocalization column.
"""

import numpy as np
import pandas as pd
from scipy.interpolate import interpn
from zarrnii import ZarrNii

# Load the regionprops table
df = pd.read_parquet(snakemake.input.regionprops_parquet)

coord_cols = snakemake.params.coord_column_names
mask_stain = snakemake.params.mask_stain
zarrnii_kwargs = snakemake.params.zarrnii_kwargs

# Validate coordinate columns
for col in coord_cols:
    if col not in df.columns:
        raise ValueError(
            f"Coordinate column '{col}' not found in regionprops. "
            f"Available columns: {df.columns.tolist()}"
        )
if len(coord_cols) != 3:
    raise ValueError(
        f"Expected exactly 3 coordinate columns, got {len(coord_cols)}: {coord_cols}"
    )

# Validate that nvoxels column is present for radius calculation
if "nvoxels" not in df.columns:
    raise ValueError(
        f"Column 'nvoxels' not found in regionprops. "
        f"Available columns: {df.columns.tolist()}"
    )

# Construct the output column name
coloc_col = f"sdt_{mask_stain}"

# Handle empty regionprops
if len(df) == 0:
    df[coloc_col] = pd.Series(dtype=float)
    df.to_parquet(snakemake.output.maskcoloc_parquet, index=False)
    print(f"Empty regionprops – wrote 0-row parquet with column '{coloc_col}'.")
else:
    # Load the signed distance transform at the finest available level (level=0).
    # The SDT voxel values are in physical units (matching the OME-Zarr coordinate
    # transformations), so the affine matrix maps voxel indices to physical space.
    znimg = ZarrNii.from_ome_zarr(snakemake.input.dist, level=0, **zarrnii_kwargs)

    # Compute physical voxel volume for the radius calculation.
    # znimg.scale is a dict keyed by dimension name, e.g. {'z': 0.004, 'y': 0.0027, 'x': 0.0027}
    scale = znimg.scale
    spatial_dims = ["z", "y", "x"]
    spacing = np.array([scale.get(d, 1.0) for d in spatial_dims], dtype=float)
    voxel_volume = float(np.prod(spacing))

    # Physical centroid coordinates from the regionprops table (N × 3, XYZ order)
    coords_xyz = df[coord_cols].values.astype(float)
    n_points = len(coords_xyz)

    # -----------------------------------------------------------------------
    # Sample SDT values at voxel coordinates via linear interpolation
    # -----------------------------------------------------------------------
    sdt_values = znimg.sample_at_points(coords_xyz)

    # -----------------------------------------------------------------------
    # Compute equivalent sphere radius from nvoxels.
    # Sphere volume:  V = (4/3) π r³
    # Solving for r:  r = ∛( 3V / (4π) )
    # -----------------------------------------------------------------------
    _SPHERE_COEFF = 3.0 / (4.0 * np.pi)  # coefficient in r = ∛(3V / (4π))
    nvoxels = df["nvoxels"].values.astype(float)
    radius = np.cbrt(_SPHERE_COEFF * nvoxels * voxel_volume)

    # -----------------------------------------------------------------------
    # Colocalization score: sdt_value - radius
    #   ≈ 0  : instance sphere boundary just touches the mask boundary
    #   < 0  : instance overlaps with the mask
    #   > 0  : instance is completely outside and away from the mask
    # -----------------------------------------------------------------------
    df[coloc_col] = sdt_values - radius

    df.to_parquet(snakemake.output.maskcoloc_parquet, index=False)

    n_overlap = int((df[coloc_col] < 0).sum())
    n_nan = int(df[coloc_col].isna().sum())
    print(f"Mask colocalization ({mask_stain}) complete:")
    print(f"  Instances analysed : {n_points}")
    print(f"  Overlapping (< 0)  : {n_overlap}")
    print(f"  Out-of-bounds (NaN): {n_nan}")
    print(f"  Output column      : '{coloc_col}'")
