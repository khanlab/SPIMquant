from zarrnii import ZarrNii
from dask.diagnostics import ProgressBar
from lib.utils import get_zarr_store, get_channel_index


store = get_zarr_store(snakemake.params.uri)
channel_index = get_channel_index(store, snakemake.wildcards.stain)

znimg = ZarrNii.from_ome_zarr(
    store, level=int(snakemake.wildcards.level), channels=[channel_index]
)


# before updating zarrnii ngffzarr3 branch to accommodate anisotropically downsampled data, instead
# we will calculate the z downsampling factor and downsample accordingly - TODO: move this to zarrnii

import numpy as np

# Get scale and axes order
scale = znimg.coordinate_transformations[0].scale

axes = znimg.axes  # list of Axis objects

# Build a mapping from axis name to index
axis_index = {axis.name.lower(): i for i, axis in enumerate(axes)}

# Extract x and z scales
x_scale = scale[axis_index["x"]]
z_scale = scale[axis_index["z"]]

# Compute ratio and power
ratio = x_scale / z_scale
level = int(np.log2(round(ratio)))


with ProgressBar():
    if level == 0:
        znimg.to_nifti(snakemake.output.nii)
    else:
        znimg.downsample(along_z=2**level).to_nifti(snakemake.output.nii)
