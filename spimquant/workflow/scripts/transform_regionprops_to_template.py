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

# Extract the coordinates as a numpy array (N x 3) in physical RAS coordinates
points = df[coord_cols].values

# Load the RAS affine transform from greedy registration
# According to resample_labels_to_zarr.py: "affine_inv_xfm = np.linalg.inv(np.loadtxt(in_xfm))"
# is used to go from spim ras to template ras. So the xfm_ras file is template->subject.
# For transforming points from subject to template, we invert it.
xfm_ras = np.loadtxt(snakemake.input.xfm_ras)
affine_subj_to_template = np.linalg.inv(xfm_ras)

# Apply affine transform to points (add homogeneous coordinate)
points_homog = np.hstack([points, np.ones((points.shape[0], 1))])  # N x 4
points_after_affine = (affine_subj_to_template @ points_homog.T).T[:, :3]  # N x 3

# Load the inverse warp (template-to-subject warp)
# When transforming points (not images), we use the opposite warp direction
# The invwarp is template->subject displacement field
warp_nib = nib.load(snakemake.input.invwarp)
warp_data = warp_nib.get_fdata()  # 4D array (X, Y, Z, 3) with displacement vectors
warp_affine = warp_nib.affine

# Convert affine-transformed points to warp's voxel coordinates for interpolation
warp_ras2vox = np.linalg.inv(warp_affine)
points_warp_vox_homog = np.hstack([points_after_affine, np.ones((points_after_affine.shape[0], 1))])
points_warp_vox = (warp_ras2vox @ points_warp_vox_homog.T).T[:, :3]

# Interpolate the warp displacement at these voxel coordinates
# map_coordinates expects coordinates in array index order
displacements = np.zeros_like(points_warp_vox)
for i in range(3):
    # Interpolate the i-th component of the displacement field
    # Note: warp_data has shape (X, Y, Z, 3) and we need (i, j, k) indexing
    coords = points_warp_vox.T  # Transpose to (3, N) for map_coordinates
    displacements[:, i] = map_coordinates(
        warp_data[..., i], coords, order=1, mode='constant', cval=0
    )

# The greedy warp is a displacement field in physical (RAS) coordinates
# So we add the displacement directly to the RAS coordinates
points_transformed = points_after_affine + displacements

# Add the transformed coordinates as new columns
df['template_x'] = points_transformed[:, 0]
df['template_y'] = points_transformed[:, 1]
df['template_z'] = points_transformed[:, 2]

# Save the output
df.to_parquet(snakemake.output.regionprops_transformed_parquet, index=False)
