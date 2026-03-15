"""Sample scalars at points, e.g. to get signed distance transform

Input (via snakemake.input):
    parquet: Parquet file with per-instance regionprops including
        centroid columns (pos_x, pos_y, pos_z).
    scalar: OME-Zarr directory containing the scalar map to sample

Parameters (via snakemake.params):
    coord_column_names: List of three column names for the x, y, z centroid
        coordinates in the regionprops table (e.g. ['pos_x', 'pos_y', 'pos_z']).
    col_name: Name of the new column
    zarrnii_kwargs: Additional keyword arguments passed to ZarrNii.from_ome_zarr
        (e.g. {'orientation': 'RPI'}).

Output (via snakemake.output):
    parquet: Parquet file identical to the input parquet but with
        the additional new column with sampled scalar
"""

import numpy as np
import pandas as pd
from zarrnii import ZarrNii

# Load the table
df = pd.read_parquet(snakemake.input.parquet)

coord_cols = snakemake.params.coord_column_names
col_name = snakemake.params.col_name
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


# Handle empty regionprops
if len(df) == 0:
    df[col_name] = pd.Series(dtype=float)
    df.to_parquet(snakemake.output.parquet, index=False)
    print(f"Empty regionprops – wrote 0-row parquet with column '{col_name}'.")
else:
    znimg = ZarrNii.from_ome_zarr(snakemake.input.scalar, **zarrnii_kwargs)

    # Physical centroid coordinates from the regionprops table (N × 3, XYZ order)
    coords_xyz = df[coord_cols].values.astype(float)
    n_points = len(coords_xyz)

    # -----------------------------------------------------------------------
    # Sample scalar values at voxel coordinates via linear interpolation
    # -----------------------------------------------------------------------
    values = znimg.sample_at_points(coords_xyz).reshape((n_points))

    df[col_name] = values

    df.to_parquet(snakemake.output.parquet, index=False)
