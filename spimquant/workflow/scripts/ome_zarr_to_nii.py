import zarr
import json
import numpy as np
import nibabel as nib
import dask.array as da
from math import log2, floor
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
from dask.diagnostics import ProgressBar


in_zarr = snakemake.input.zarr

zi = zarr.open(in_zarr)

attrs=zi['/'].attrs.asdict()

#get channel index from omero metadata
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)

level=int(snakemake.wildcards.level)

#read coordinate transform from ome-zarr
transforms = attrs['multiscales'][0]['datasets'][level]['coordinateTransformations']


#zarr uses z,y,x ordering, we reverse this for nifti
# also flip to set orientation properly
affine = np.eye(4)
affine[0,0]=-transforms[0]['scale'][-1] #x
affine[1,1]=-transforms[0]['scale'][-2] #y
affine[2,2]=-transforms[0]['scale'][-3] #z


#grab the channel index corresponding to the stain
darr = da.from_zarr(in_zarr,component=f'/{level}')[channel_index,:,:,:].squeeze()

#downsample in z
#calculate scaling between x and z -- if scaling is higher in z, then leave as is.. 
#if z is smaller, then we downsample by power of 2, one less than what would make it greater than x
z_ratio = transforms[0]['scale'][-1] / transforms[0]['scale'][-3]  # x / z scaling (if >1, z needs to be downsampled)
zdownsampling = 2**(floor(log2(z_ratio)))

if zdownsampling > 1:

    affine[2,2]=zdownsampling*affine[2,2]

    #  we achieve this by rechunking, then performing mean over axis-0 in each block
    in_chunksize=(zdownsampling,darr.shape[1],darr.shape[2])
    out_chunksize=(1,darr.shape[1],darr.shape[2])

    darr = darr.rechunk(in_chunksize).map_blocks(lambda x: np.mean(x,axis=0).reshape(1,x.shape[1],x.shape[2]),chunks=out_chunksize)


#input array axes are ZYX 
#writing to nifti we want XYZ
darr_mvax = da.moveaxis(darr,(0,1,2),(2,1,0))

with ProgressBar():
    nii = nib.Nifti1Image(darr_mvax.compute(),
                    affine=affine
                    )
                    
nii.to_filename(snakemake.output.nii)
