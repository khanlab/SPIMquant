import zarr
import dask.array as da
from dask.array.overlap import overlap
import numpy as np
from math import sqrt
from skimage.feature import blob_dog
from dask.diagnostics import ProgressBar


in_zarr = snakemake.input.zarr

zi = zarr.open(in_zarr)
attrs=zi['/'].attrs.asdict()

#get channel index from omero metadata
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)


level=snakemake.params.level

#read coordinate transform from ome-zarr
transforms = attrs['multiscales'][0]['datasets'][level]['coordinateTransformations']


#darr_chan = da.from_zarr(in_zarr,component=f'{level}',chunks=snakemake.params.chunks)[channel_index,:,:,:]
darr_chan = da.from_zarr(in_zarr,component=f'{level}',chunks=(1,1,50,50))[channel_index,:,:,:]

print(f'shape: {darr_chan.shape}')
print(f'chunks: {darr_chan.chunks}')

#adjust sigma based on physical size of voxels --- TODO check this! 
# get um per pixel from the scaling transforms
# then convert from um to pixel by dividing by it
scaling_zyx=np.array(transforms[0]['scale'][1:])
print(f'scaling_zyx: {scaling_zyx}') 
#mm per pixel

min_sigma_px = snakemake.params.min_sigma_um *1e-3 / scaling_zyx 
max_sigma_px = snakemake.params.max_sigma_um *1e-3 / scaling_zyx
boundary_px = tuple(max_sigma_px.astype('int').tolist())

print(f'min_sigma_px: {min_sigma_px}')
print(f'max_sigma_px: {max_sigma_px}')
print(f'boundary_px: {boundary_px}')


def detect_blobs(x,block_info=None):

    #we need local chunk location in order to translate the local blob
    #coords into global blob coords
    arr_location = block_info[0]['array-location']

    blobs_dog = blob_dog(x, min_sigma=min_sigma_px, max_sigma=max_sigma_px, 
            threshold=snakemake.params.threshold, exclude_border=boundary_px)
    blobs_dog[:, -1] = blobs_dog[:, -1] * (2 ** 0.3333)  #adjust to get radius

    #offset by chunk location, then scale by header
    for ax in range(3):
        #scale the coordinates
        blobs_dog[:,ax] = scaling_zyx[ax] * (arr_location[ax][0] + blobs_dog[:,ax])

        #and the radii 
        blobs_dog[:,ax+3] =  scaling_zyx[ax] * blobs_dog[:,ax+3]


    return blobs_dog

#expanded = overlap(darr_chan, depth=boundary_px, boundary=0)

#darr_blobs = expanded.map_blocks(detect_blobs,drop_axis=[2],dtype='float')
darr_blobs = darr_chan.map_blocks(detect_blobs,drop_axis=[2],dtype='float')

with ProgressBar():
#    da.to_zarr(darr_blobs,snakemake.output.zarr)
    computed_blobs = darr_blobs.compute()

#TODO: pick a better format? tsv? 
np.save(snakemake.output.npy,computed_blobs)

