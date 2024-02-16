import zarr
import dask.array as da
import nibabel as nib
import numpy as np
from scipy.interpolate import interpn
from dask.diagnostics import ProgressBar



in_zarr=snakemake.input.zarr_zip
in_xfm=snakemake.input.xfm_ras
in_template_dseg=snakemake.input.dseg
out_zarr=snakemake.output.zarr

level=0

#load dask array from zarr reference image
darr = da.from_zarr(in_zarr,component=f'/{level}')

#load template dseg
dseg_nib = nib.load(in_template_dseg)

# get dseg volume for interpolation,
dseg_vol = dseg_nib.get_fdata()

# along with grid points for interpolation
grid_points = (np.arange(dseg_vol.shape[0]),
             np.arange(dseg_vol.shape[1]),
             np.arange(dseg_vol.shape[2]))



#read coordinate transform from ome-zarr
zi = zarr.open(in_zarr)
attrs=zi['/'].attrs.asdict()
level=0
transforms = attrs['multiscales'][0]['datasets'][level]['coordinateTransformations']

#need to put together the sequence of transforms to apply

# 1. scaling_xfm (vox2ras in spim space)
#note: zarr uses z,y,x ordering, we reverse this for nifti; also the negation to flip
# this matches what the ome_zarr_to_nii affine has
scaling_xfm = np.eye(4)
scaling_xfm[0,0]=-transforms[0]['scale'][2] #x 
scaling_xfm[1,1]=-transforms[0]['scale'][1] #y
scaling_xfm[2,2]=-transforms[0]['scale'][0] #z

# 2. affine_inv_xfm (from registration, takes points from spim ras to template ras)
affine_inv_xfm = np.linalg.inv(np.loadtxt(in_xfm))


# 3. ras2vox in template space
ras2vox = np.linalg.inv(dseg_nib.affine)

# concatenate all three
concat_xfm = ras2vox @ affine_inv_xfm @ scaling_xfm



def interp_label(x,block_info=None):
   # print(block_info)

    arr_location = block_info[0]['array-location']
    
    xv,yv,zv=np.meshgrid(np.arange(arr_location[2][0],arr_location[2][1]),
            np.arange(arr_location[1][0],arr_location[1][1]),
            np.arange(arr_location[0][0],arr_location[0][1]))

    #print(f'x shape: {x.shape}')
    #print(f'xv shape: {xv.shape}')

    #reshape them into a vectors (x,y,z,1) for each point, so we can matrix multiply
    xvf=xv.reshape((1,np.product(xv.shape)))
    yvf=yv.reshape((1,np.product(yv.shape)))
    zvf=zv.reshape((1,np.product(zv.shape)))
    homog=np.ones(xvf.shape)
    
    vecs=np.vstack((xvf,yvf,zvf,homog))
    
    xfm_vecs = concat_xfm @ vecs
    
    #then finally interpolate those points on the template dseg volume
    interpolated = interpn(grid_points,dseg_vol,
                        xfm_vecs[:3,:].T, #
                        method='nearest',
                        bounds_error=False,
                        fill_value=-100)
    
    return interpolated.reshape(x.shape)


#perform interpolation on each block of spim zarr, in parallel
darr_map=darr.map_blocks(interp_label, dtype=np.uint16)

with ProgressBar():
    da.to_zarr(darr_map,out_zarr,overwrite=True)


