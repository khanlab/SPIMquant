from scipy.interpolate import interpn
import nibabel as nib
import numpy as np

# we have coordinates of blobs, along with an image with the boundary distance penalty

# we just need to perform an interpolation of each point into the volume, and include/exclude based on a threshold

# in the future we could make use of additional features (like eigenvalues of hessian) as well 


# load the nii brainmask penalty
penalty_nib = nib.load(snakemake.input.penalty)
penalty_vol = penalty_nib.get_fdata()
penalty_ras2vox = np.linalg.inv(penalty_nib.affine) 

print(f'ras2vox: {penalty_ras2vox}')
print(f'vox2ras: {penalty_nib.affine}')

penalty_grid = (np.arange(penalty_vol.shape[0]),
        np.arange(penalty_vol.shape[1]),
        np.arange(penalty_vol.shape[2]))

# load the coordinates of points (Nx3)
points = np.load(snakemake.input.points_npy)

print(f'points, mean: {points.mean(axis=0)}, max: {points.max(axis=0)}, min: {points.min(axis=0)}')

# points are in RAS (physical) space -- we want to transform them to voxel 
# space (relative to the penalty_vol), so we can use the penalty_ras2vox for this

# but we flip (negate each dim) when making nifti, so this seemingly needs to be done to the points before transforming
# we definitely still need to reverse the ordering from z y x to x y z
# -- TODO look into this..
#  one way to keep it internally consistent, is to make sure we use ome-zarr-neuro to get the scalings always --
#  e.g. get_vox2ras_zarr (which includes reordering and flipping)
#    that code will need to be revised to make sure the ome-zarr display space is also included (ie for napari visualization), which doesn't have any reordering or flipping

points_ras = -np.flip(points,axis=1)

points_vox = penalty_ras2vox @ np.hstack((points_ras,np.ones((points_ras.shape[0],1)))).T

print(f'points_vox, mean: {points_vox.mean(axis=1)}, max: {points_vox.max(axis=1)}, min: {points_vox.min(axis=1)}')

# now that we have points in voxel space, we can use interpn to interpolate
penalty_at_points = interpn(penalty_grid,
                                penalty_vol,
                                points_vox[:3,:].T,
                                method='linear',
                                bounds_error=False,
                                fill_value=0)

print(f'penalty at points, mean: {penalty_at_points.mean()}, max: {penalty_at_points.max()}, min: {penalty_at_points.min()}')
#now I have the Nx3 points, and a Nx1 penalty_at_points

#save the filtered points
filtered_points = points[penalty_at_points > snakemake.params.threshold,:]

n_before = points.shape[0]
n_after= filtered_points.shape[0]

percent_removed = 100.0* (n_before - n_after)/n_before
print(f'filtering removed {percent_removed:02f}% of points, from {n_before} to {n_after}')

np.save(snakemake.output.points_npy,filtered_points)


