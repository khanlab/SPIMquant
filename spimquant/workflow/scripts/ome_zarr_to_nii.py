import zarr
import json
import numpy as np
import nibabel as nib
import dask.array as da
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader

#TODO: grab channel name from OMERO metadata 

in_zarr = snakemake.input.zarr
channel_index = snakemake.params.channel_index

zi = zarr.open(in_zarr)

attrs=zi['/'].attrs.asdict()

level=int(snakemake.wildcards.level)

#read coordinate transform from ome-zarr
transforms = attrs['multiscales'][0]['datasets'][level]['coordinateTransformations']


#zarr uses z,y,x ordering, we reverse this for nifti
# also flip to set orientation properly
affine = np.eye(4)
affine[0,0]=-transforms[0]['scale'][3] #x
affine[1,1]=-transforms[0]['scale'][2] #y
affine[2,2]=-transforms[0]['scale'][1] #z

#grab the channel index corresponding to the stain
darr = da.from_zarr(in_zarr,component=f'/{level}')[channel_index,:,:,:].squeeze()

#input array axes are ZYX 
#writing to nifti we want XYZ
out_arr = np.moveaxis(darr,(0,1,2),(2,1,0)) 

nii = nib.Nifti1Image(out_arr,
                    affine=affine
                    )
                    
nii.to_filename(snakemake.output.nii)
