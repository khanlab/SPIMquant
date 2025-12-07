"""Transform regionprops coordinates from subject to template space.

This script takes computed regionprops from segmentation in subject space
and applies the affine and warp transforms to convert the spatial coordinates
(pos_x, pos_y, pos_z) from subject space to template space, adding new columns
(template_x, template_y, template_z) while retaining all other regionprops columns.
"""

import pandas as pd
import numpy as np
import nibabel as nib
from scipy.ndimage import map_coordinates

# Load the regionprops data
df = pd.read_parquet(snakemake.input.regionprops_parquet)

# Get coordinate column names from config
coord_cols = snakemake.params.coord_column_names  # e.g., ['pos_x', 'pos_y', 'pos_z']

# Validate that coordinate columns exist
for col in coord_cols:
    if col not in df.columns:
        raise ValueError(
            f"Coordinate column '{col}' not found in regionprops data. Available columns: {df.columns.tolist()}"
        )

# Extract the coordinates as a numpy array (N x 3) in physical RAS coordinates
points = df[coord_cols].values

# Validate that coordinates are numeric
if not np.issubdtype(points.dtype, np.number):
    raise ValueError(
        f"Coordinate columns must contain numeric data, got dtype: {points.dtype}"
    )


# Transform points by applying inverse affine, then inverse warp
invaff = AffineTransform.from_txt(snakemake.input.xfm_ras).invert()
invwarp = DisplacementTransform.from_nifti(snakemake.input.invwarp)


points_transformed = invwarp.apply_transform(invaff.apply_transform(points.T))


# Add the transformed coordinates as new columns
df["template_x"] = points_transformed[:, 0]
df["template_y"] = points_transformed[:, 1]
df["template_z"] = points_transformed[:, 2]

# Save the output
df.to_parquet(snakemake.output.regionprops_transformed_parquet, index=False)
