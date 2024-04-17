import dask.array as da
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


darr_chan = da.from_zarr(in_zarr,component=f'{level}')[channel_index,:,:,:]

#adjust sigma based on physical size of voxels --- TODO check this! 
raw_pix_um=transforms[0]['scale'][1:]

ds_byax = np.array((1,2**level,2**level))
ds_um_per_pix = raw_pix_um*ds_byax
ds_um_per_pix

min_sigma_px = min_sigma_um / ds_um_per_pix
max_sigma_px = max_sigma_um / ds_um_per_pix




def detect_blobs(x,block_info=None):

    #we need local chunk location in order to translate the local blob
    #coords into global blob coords
    arr_location = block_info[0]['array-location']

    blobs_dog = blob_dog(x, min_sigma=min_sigma_px, max_sigma=max_sigma_px, threshold=threshold)
    blobs_dog[:, -1] = blobs_dog[:, -1] * sqrt(2) #adjust to get radius

    #offset by chunk location, then scale by header
    for ax in range(3):
        #scale the coordinates

        #TODO: add exclude on edges?
        blobs_dog[:,ax] = scaling_zyx[ax] * (arr_location[ax][0] + blobs_dog[:,ax])

        #and the radii #TODO CHECK THIS!
        blobs_dog[:,ax+3] =  (arr_location[ax][0] + blobs_dog[:,ax+3]) / scaling_zyx[ax]


    return blobs_dog

#TODO: pick depth according to max sigma? 
expanded = overlap(arr, depth=2, boundary=0)

darr_blobs = darr_chan.map_blocks(detect_blobs,drop_axis=[2],dtype='float')

with ProgressBar():
    computed_blobs = darr_blobs.compute()

#TODO: pick a better format? napari points layer native?
np.save(snakemake.output.npy,computed_blobs)

