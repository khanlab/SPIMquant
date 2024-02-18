import zarr
import dask.array as da
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.interpolate import interpn
from dask.diagnostics import ProgressBar
from ome_zarr.writer import write_labels, write_label_metadata
from ome_zarr.scale import Scaler




in_zarr=snakemake.input.zarr_zip
in_xfm=snakemake.input.xfm_ras
in_template_dseg=snakemake.input.dseg
out_zarr=snakemake.output.zarr
level=snakemake.params.level_to_resample_to
max_layer=snakemake.params.max_downsampling_layers #number of downsamplings by 2 to include in zarr
label_tsv = snakemake.input.label_tsv
scaling_method=snakemake.params.scaling_method
label_name=snakemake.params.label_name



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
multiscale=0 #first multiscale image
transforms = attrs['multiscales'][multiscale]['datasets'][level]['coordinateTransformations']

#for writing metadata:
axes = attrs['multiscales'][multiscale]['axes']
coordinate_transformations = [ attrs['multiscales'][multiscale]['datasets'][level]['coordinateTransformations'] for level in range(level,max_layer+1)]


#need to put together the sequence of transforms to apply

# 1. scaling_xfm (vox2ras in spim space)
#note: zarr uses z,y,x ordering, we reverse this for nifti; also the negation to flip
# this matches what the ome_zarr_to_nii affine has
scaling_xfm = np.eye(4)
scaling_xfm[0,0]=-transforms[0]['scale'][2] #x  # 0-index in transforms is the first (and only) transform 
scaling_xfm[1,1]=-transforms[0]['scale'][1] #y
scaling_xfm[2,2]=-transforms[0]['scale'][0] #z

# 2. affine_inv_xfm (from registration, takes points from spim ras to template ras)
affine_inv_xfm = np.linalg.inv(np.loadtxt(in_xfm))


# 3. ras2vox in template space
ras2vox = np.linalg.inv(dseg_nib.affine)

# concatenate all three
concat_xfm = ras2vox @ affine_inv_xfm @ scaling_xfm



def interp_label(x,block_info=None):

    arr_location = block_info[0]['array-location']
    
    xv,yv,zv=np.meshgrid(np.arange(arr_location[2][0],arr_location[2][1]),
            np.arange(arr_location[1][0],arr_location[1][1]),
            np.arange(arr_location[0][0],arr_location[0][1]))


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
                        fill_value=0)
    
    return interpolated.reshape(x.shape)


#perform interpolation on each block of spim zarr, in parallel
darr_map=darr.map_blocks(interp_label, dtype=np.uint16)

#write ome-zarr metadata
store = zarr.DirectoryStore(out_zarr,dimension_separator='/')

#first, copy entire input image (ie channel data) to the output zarr, 
# as that is needed alongside the labels to visualize
root_in = zi['/']
root = zarr.group(store,path='/',overwrite=True)

#then we add in labels
scaler = Scaler(max_layer=max_layer-level,method='nearest')

#label metadata - convert from bids tsv to list of dicts
df = pd.read_csv(label_tsv,sep='\t',index_col='index')

colors=[]
properties=[]
for i, row in df.iterrows():
    if i>0: #ensure bg label is transparent
        alpha=255
    else:
        alpha=0

    hexrgb=row['color'].lstrip("#")
    rgba = [int(hexrgb[i:i+2], 16) for i in (0, 2, 4)]
    rgba.append(alpha)
    colors.append({'label-value': i, 'rgba': tuple(rgba)})
    properties.append({'label-value': i, 'name': row['name'], 'abbreviation': row['abbreviation']})  

with ProgressBar():
    zarr.copy_all(root_in,root)
    write_labels(labels=darr_map,
                            group=root,
                            scaler=scaler,
                            name=label_name,
                            coordinate_transformations=coordinate_transformations,
                            storage_options={'dimension_separator': '/'},
                            axes=axes)

    write_label_metadata(group=root,
                            name=f'/labels/{label_name}',
                            colors=colors,
                            properties=properties)

