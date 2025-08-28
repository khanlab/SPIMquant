from zarrnii import ZarrNii
from lib.utils import get_zarr_store, get_channel_index


store = get_zarr_store(snakemake.params.spim_uri)

channel_index = get_channel_index(store, snakemake.wildcards.stain)

in_orient = snakemake.config[
    "in_orientation"
]  # TODO: this is a bit ugly - update the ZarrNii to recognize None  as an orientation, so we can just pass in_orientation
orient_opt = {} if in_orient == None else {"orientation": in_orient}

hires_level = int(snakemake.wildcards.level)
ds_level = int(snakemake.wildcards.dslevel)

# this function optionally uses coiled


znimg_hires = ZarrNii.from_ome_zarr(
    store, channels=[channel_index], level=hires_level, **orient_opt
)

# --------------
# before updating zarrnii ngffzarr3 branch to accommodate anisotropically downsampled data, instead
# we will calculate the z downsampling factor and downsample accordingly - TODO: move this to zarrnii

import numpy as np

# Get scale and axes order
scale = znimg_hires.coordinate_transformations[0].scale

axes = znimg_hires.axes  # list of Axis objects

# Build a mapping from axis name to index
axis_index = {axis.name.lower(): i for i, axis in enumerate(axes)}

# Extract x and z scales
x_scale = scale[axis_index["x"]]
z_scale = scale[axis_index["z"]]

# Compute ratio and power
ratio = x_scale / z_scale
level = int(np.log2(round(ratio)))

# ----------


if level > 0:
    znimg_hires = znimg_hires.downsample(along_z=2**level)

hires_shape = znimg_hires.darr.shape


# ok, now we have the bias field.. let's resample it to the level
# where we perform thresholding

# Compute step 1: first, we write the nifti n4 bias field to ome zarr (this is so we can use distributed computing)
ZarrNii.from_nifti(snakemake.input.n4_bf_ds, chunks=(1, 10, 10, 10)).to_ome_zarr(
    snakemake.params.bf_ds_uri
)

if snakemake.config["use_coiled"]:
    from coiled import Cluster

    cluster = Cluster(
        name="coiled-snakemake", package_sync_ignore=["spimquant"], n_workers=10
    )
    client = cluster.get_client()


# Compute step 2: then we upsample the bias field, and save to file
znimg_biasfield_upsampled = ZarrNii.from_ome_zarr(
    snakemake.params.bf_ds_uri, **orient_opt
).upsample(to_shape=hires_shape)
znimg_biasfield_upsampled.to_ome_zarr(snakemake.params.bf_us_uri, max_layer=0)


# Compute step 3: now multiply biasfield and hires image together.  we read the biasfield
# in lazily again (this is so the subsequent step is a different dask computation than the previous step),
# and write the full-res N4-corrected OME-zarr dataset
znimg_biasfield_upsampled = ZarrNii.from_ome_zarr(
    snakemake.params.bf_us_uri, **orient_opt
)

if level == 0:

    znimg_hires = ZarrNii.from_ome_zarr(
        store,
        channels=[channel_index],
        level=hires_level,
        chunks=znimg_biasfield_upsampled.darr.chunks,
        rechunk=True,
        **orient_opt,
    )  # chunk size comes from upsampled array

else:
    znimg_hires = ZarrNii.from_ome_zarr(
        store,
        channels=[channel_index],
        level=hires_level,
        **orient_opt,
    ).downsample(along_z=2**level)
    znimg_hires.darr = znimg_hires.darr.rechunk(znimg_biasfield_upsampled.darr.chunks)


znimg_hires.darr = znimg_hires.darr / znimg_biasfield_upsampled.darr
znimg_hires.to_ome_zarr(snakemake.params.spim_n4_uri, max_layer=5)
