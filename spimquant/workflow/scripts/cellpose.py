import relabel
import numpy as np
import zarr
import dask.array as da
from dask.diagnostics import ProgressBar
from skimage import color
import torch
from cellpose import models


# ## Load the 3D Cellpose model

model = models.Cellpose(gpu=False, model_type='nuclei')


# see https://cellpose.readthedocs.io/en/latest/settings.html

#diameter is an important parameter - set to expected size of cells
def cellposeLabel(im_chunk, model=None):
    masks, _, _, _ = model.eval(im_chunk, diameter=64, channels=[0, 0], flow_threshold=0.4, do_3D=True)
    masks = masks.astype(np.int32)
    return masks

in_zarr = snakemake.input.zarr
zi = zarr.open(in_zarr)
attrs=zi['/'].attrs.asdict()

#get channel index from omero metadata
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)

level = snakemake.params.level

img_da = da.from_zarr(in_zarr,component=f'{level}',chunks=snakemake.params.chunks)[channel_index,:,:,:]


labels = relabel.image2labels(
        img_da,seg_fn=cellposeLabel,
        overlaps=[0, 64, 64],
        ndim=3,
        segmentation_fn_kwargs={"model": model})


with ProgressBar():
    labels.rechunk(chunks=snakemake.params.chunks).to_zarr(snakemake.output.zarr)




