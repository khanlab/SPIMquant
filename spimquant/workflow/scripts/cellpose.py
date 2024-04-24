import relabel
import numpy as np
import zarr
import dask.array as da
from dask.diagnostics import ProgressBar
from skimage import color
import torch
from cellpose import models
from dask.distributed import LocalCluster
import dask

dask.config.set(scheduler='threads',num_workers=snakemake.threads)


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

img_da = da.from_zarr(in_zarr,component=f'{level}')[channel_index,:,:,:].squeeze().rechunk(chunks=snakemake.params.chunks)

labels = da.map_blocks(
        cellposeLabel,
        img_da,
        model=model,
        dtype=np.int32,
        meta=np.empty((0, 0), dtype=np.int32)
    )


with ProgressBar():
    labels.to_zarr(snakemake.output.zarr)

