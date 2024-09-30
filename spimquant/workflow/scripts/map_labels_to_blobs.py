from scipy.interpolate import interpn
import nibabel as nib
import numpy as np
import pandas as pd

# load the nii dseg
dseg_nib = nib.load(snakemake.input.dseg)
dseg_vol = dseg_nib.get_fdata()
dseg_ras2vox = np.linalg.inv(dseg_nib.affine) 

print(f'ras2vox: {dseg_ras2vox}')
print(f'vox2ras: {dseg_nib.affine}')

dseg_grid = (np.arange(dseg_vol.shape[0]),
        np.arange(dseg_vol.shape[1]),
        np.arange(dseg_vol.shape[2]))

# load the coordinates of points (Nx3)
points = np.load(snakemake.input.points_npy)

if points.shape[0] == 0:
    print('No points found, skip points statistics')
else:
    print(f'points, shape: {points.shape}, mean: {points.mean(axis=0)}, max: {points.max(axis=0)}, min: {points.min(axis=0)}')

# points are in RAS (physical) space -- we want to transform them to voxel 
# space (relative to the dseg_vol), so we can use the dseg_ras2vox for this

# but we flip (negate each dim) when making nifti, so this seemingly needs to be done to the points before transforming
# we definitely still need to reverse the ordering from z y x to x y z
# -- TODO look into this..
#  one way to keep it internally consistent, is to make sure we use ome-zarr-neuro to get the scalings always --
#  e.g. get_vox2ras_zarr (which includes reordering and flipping)
#    that code will need to be revised to make sure the ome-zarr display space is also included (ie for napari visualization), which doesn't have any reordering or flipping

points_ras = -np.flip(points,axis=1)

points_vox = dseg_ras2vox @ np.hstack((points_ras,np.ones((points_ras.shape[0],1)))).T

if points.shape[0] != 0:
    print(f'points_vox, mean: {points_vox.mean(axis=1)}, max: {points_vox.max(axis=1)}, min: {points_vox.min(axis=1)}')

# now that we have points in voxel space, we can use interpn to interpolate
label_at_points = interpn(dseg_grid,
                                dseg_vol,
                                points_vox[:3,:].T,
                                method='nearest',
                                bounds_error=False,
                                fill_value=0)

if points.shape[0] != 0:
    print(f'label_at_points, shape: {label_at_points.shape}, mean: {label_at_points.mean()}, max: {label_at_points.max()}, min: {label_at_points.min()}')

#now I have the Nx3 points, and a Nx1 label_at_points


# TODO instead of concatenating arrays, should concat dataframes (to keep points as float, and labels as int)

# TODO add label names too (lookup from dseg tsv)


#save the points plus labels as a tsv
df = pd.DataFrame(np.hstack((points,label_at_points.reshape((label_at_points.shape[0],1)))),columns=['z','y','x','label_index'])



df.to_csv(snakemake.output.blobs_tsv,sep='\t',index=False)



