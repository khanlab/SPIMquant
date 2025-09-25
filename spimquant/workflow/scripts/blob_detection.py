from math import sqrt

import dask
import dask.array as da
import numpy as np
import sparse
import zarr
from dask.array.overlap import overlap, trim_overlap
from dask.diagnostics import ProgressBar
from skimage.feature import blob_dog

# set threads
dask.config.set(scheduler="threads", num_workers=snakemake.threads)

in_zarr = snakemake.params.zarr

zi = zarr.open(in_zarr)
attrs = zi["/"].attrs.asdict()

# get channel index from omero metadata
channel_labels = [channel_dict["label"] for channel_dict in attrs["omero"]["channels"]]
channel_index = channel_labels.index(snakemake.wildcards.stain)


level = snakemake.params.level

# read coordinate transform from ome-zarr
transforms = attrs["multiscales"][0]["datasets"][level]["coordinateTransformations"]


darr_chan = da.from_zarr(in_zarr, component=f"{level}", chunks=snakemake.params.chunks)[
    channel_index, :, :, :
]

print(f"before adding overlap:")
print(f"shape: {darr_chan.shape}")
print(f"chunks: {darr_chan.chunks}")

# adjust sigma based on physical size of voxels --- TODO check this!
# get um per pixel from the scaling transforms
# then convert from um to pixel by dividing by it
scaling_zyx = np.array(transforms[0]["scale"][1:])
print(f"scaling_zyx: {scaling_zyx}")
# mm per pixel

min_sigma_px = snakemake.params.min_sigma_um * 1e-3 / scaling_zyx
max_sigma_px = snakemake.params.max_sigma_um * 1e-3 / scaling_zyx
boundary_px = tuple(max_sigma_px.astype("int").tolist())

print(f"min_sigma_px: {min_sigma_px}")
print(f"max_sigma_px: {max_sigma_px}")
print(f"boundary_px: {boundary_px}")


def detect_blobs(x, block_info=None):

    # we need local chunk location in order to translate the local blob
    # coords into global blob coords
    arr_location = block_info[0]["array-location"]

    blobs_dog = blob_dog(
        x,
        min_sigma=min_sigma_px,
        max_sigma=max_sigma_px,
        threshold=snakemake.params.threshold,
        exclude_border=boundary_px,
    )

    # we have coordinates, now convert to a volumetric representation
    # -- (index,sigma_0,sigma_1,sigma_2) -- this requires

    blobs_dog[:, 3:] = blobs_dog[:, 3:] * sqrt(
        3
    )  # adjust each sigma to get radius (this is still in pixels)

    blobs_img = np.zeros((x.shape[0], x.shape[1], x.shape[2], 4), dtype=np.float32)
    for i in range(blobs_dog.shape[0]):
        blobs_img[
            int(blobs_dog[i, 0]), int(blobs_dog[i, 1]), int(blobs_dog[i, 2]), 0
        ] = (
            i + 1
        )  # set to index for now (1-indexed)
        for sigma_i in range(3):
            blobs_img[
                int(blobs_dog[i, 0]),
                int(blobs_dog[i, 1]),
                int(blobs_dog[i, 2]),
                1 + sigma_i,
            ] = blobs_dog[i, sigma_i]

    return blobs_img


expanded = overlap(darr_chan, depth=boundary_px, boundary=0)

print(f"after adding overlap:")
print(f"shape: {expanded.shape}")
print(f"chunks: {expanded.chunks}")


darr_blobs = expanded.map_blocks(
    detect_blobs, dtype=np.float32, meta=np.array((), dtype=np.float32)
)

darr_blobs_trim = trim_overlap(darr_blobs, depth=boundary_px, boundary=0)


"""
#save to zarr 
with ProgressBar():
    darr_blobs_trim.to_zarr('temp.zarr')

#now convert the saved zarr to a sparse array
with ProgressBar():
    sparse_array = da.from_zarr('temp.zarr').map_blocks(sparse.COO).compute()
"""


# try skipping the zarr intermediary - but below might be too memory-intensive??
with ProgressBar():
    sparse_array = darr_blobs_trim.map_blocks(sparse.COO).compute()

# save the sparse array
sparse.save_npz(snakemake.output.sparse_npz, sparse_array)

# multiply coords with physical size
scaled_coords = sparse_array[:, :, :, 0].coords.T * scaling_zyx.reshape(1, 3)

np.save(snakemake.output.points_npy, scaled_coords)
