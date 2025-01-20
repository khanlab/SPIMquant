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


with ProgressBar():
    # ZarrNii now implicitly downsamples xy and z if warranted by chosen level
    ZarrNii.from_ome_zarr(store,level=level,channels=[channel_index]).to_nifti(snakemake.output.nii)

    

