import zarr
import json
import numpy as np
import nibabel as nib
import dask.array as da
from math import log2, floor
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii
from lib.cloud_io import get_fsspec, is_remote
from upath import UPath as Path

in_zarr_url = snakemake.params.in_zarr

if is_remote(in_zarr_url):
    #fs_args={'storage_provider_settings':snakemake.params.storage_provider_settings,'creds':snakemake.input.creds}
    fs_args={'creds':snakemake.input.creds}
else:
    fs_args={}

fs = get_fsspec(in_zarr_url,**fs_args)

if Path(in_zarr_url).suffix == '.zip':
    store = zarr.storage.ZipStore(Path(in_zarr_url).path,dimension_separator='/',mode='r')
else:
    store = zarr.storage.FSStore(Path(in_zarr_url).path,fs=fs,dimension_separator='/',mode='r')


zi = zarr.open(store=store,mode='r')

#zi = zarr.open(in_zarr_url)
attrs=zi['/'].attrs.asdict()

#get channel index from omero metadata
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)

level=int(snakemake.wildcards.level)


#read coordinate transform from ome-zarr
transforms = attrs['multiscales'][0]['datasets'][level]['coordinateTransformations']

x_scaling=transforms[0]['scale'][-1] #x
z_scaling=transforms[0]['scale'][-3] #z

#calculate scaling between x and z -- if scaling is higher in z, then leave as is.. 
#if z is smaller, then we downsample by power of 2, one less than what would make it greater than x
z_ratio = x_scaling / z_scaling  # x / z scaling (if >1, z needs to be downsampled)
zdownsampling = 2**(floor(log2(z_ratio)))


with ProgressBar():
    ZarrNii.from_path(store,level=level,channels=[channel_index]).downsample(along_z=zdownsampling).to_nifti(snakemake.output.nii)

    

