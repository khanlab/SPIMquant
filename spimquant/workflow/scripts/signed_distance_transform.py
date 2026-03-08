"""Compute signed distance transform from a binary mask using dask map_overlap.

The signed distance transform assigns:
- Negative values to interior voxels (inside the mask), proportional to distance
  to the nearest background voxel.
- Positive values to exterior voxels (outside the mask), proportional to distance
  to the nearest foreground voxel.

Distance values are in physical units (same units as the OME-Zarr coordinate
transformations, typically mm or µm).

Uses scipy.ndimage.distance_transform_edt via dask.array.map_overlap for chunked,
parallel processing with spatial overlap padding to reduce boundary artifacts.

For blocks that are entirely foreground or entirely background, the distance is
filled with the maximum possible distance across the block (the physical diagonal
length of the block), since the true boundary lies beyond the block extent.
"""

import dask.array as da
import numpy as np
from scipy.ndimage import distance_transform_edt
from zarrnii import ZarrNii

from dask_setup import get_dask_client

with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

    overlap_depth = snakemake.params.overlap_depth

    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.mask, **snakemake.params.zarrnii_kwargs
    )

    # Get physical voxel spacing from the ZarrNii scale metadata.
    # znimg.scale is a dict keyed by dimension name (e.g. {'z': 0.004, 'y': 0.0027, 'x': 0.0027}).
    # znimg.dims gives the ordered dimension names, e.g. ['c', 'z', 'y', 'x'].
    scale = znimg.scale
    _known_spatial = {"z", "y", "x"}
    spatial_dims = [d for d in znimg.dims if d in _known_spatial]
    spacing = []
    for d in spatial_dims:
        s = scale.get(d)
        if s is None:
            import warnings

            warnings.warn(
                f"Physical spacing for dimension '{d}' not found in OME-Zarr "
                "metadata; defaulting to 1.0. Distance values may not be in "
                "physical units.",
                stacklevel=2,
            )
            s = 1.0
        spacing.append(s)
    spacing = np.array(spacing)

    # Calculate the maximum physical distance across a single block.
    # This is used as the fill value for entirely-foreground or entirely-background
    # blocks, where the true boundary lies outside the block.
    # max() is used because the last chunk along each axis may be smaller than
    # the others (when the array size is not divisible by the chunk size); the
    # largest chunk determines the worst-case block diagonal.
    spatial_chunk_indices = [znimg.dims.index(d) for d in spatial_dims]
    block_size = np.array(
        [max(znimg.darr.chunks[i]) for i in spatial_chunk_indices]
    )
    max_dist = float(np.sqrt(np.sum((block_size * spacing) ** 2)))

    def signed_dt_block(block):
        """Compute signed distance transform for a single block.

        Expects block shape (C, Z, Y, X).  For each channel the Euclidean
        distance transform is computed twice:
          - dt_inside : distance (in physical units) from each foreground
            voxel to the nearest background voxel.
          - dt_outside: distance (in physical units) from each background
            voxel to the nearest foreground voxel.

        The signed distance transform is dt_outside - dt_inside, giving
        negative values inside the mask and positive values outside.

        Blocks that are entirely foreground are filled with -max_dist, and
        blocks that are entirely background are filled with +max_dist.
        """
        result = np.zeros(block.shape, dtype=np.float32)
        for c in range(block.shape[0]):
            binary = block[c] > 0
            n_fg = np.count_nonzero(binary)
            n_total = binary.size
            if n_fg == 0:
                # All background: nearest foreground is at least one block away
                result[c] = max_dist
            elif n_fg == n_total:
                # All foreground: nearest background is at least one block away
                result[c] = -max_dist
            else:
                dt_inside = distance_transform_edt(binary, sampling=spacing).astype(
                    np.float32
                )
                dt_outside = distance_transform_edt(
                    ~binary, sampling=spacing
                ).astype(np.float32)
                result[c] = dt_outside - dt_inside
        return result

    # depth=0 for the channel dimension, overlap_depth for spatial dims
    depth = {0: 0, 1: overlap_depth, 2: overlap_depth, 3: overlap_depth}

    sdt_darr = da.map_overlap(
        signed_dt_block,
        znimg.darr,
        depth=depth,
        boundary=0,
        dtype=np.float32,
    )

    znimg.darr = sdt_darr

    znimg.to_ome_zarr(snakemake.output.dist, max_layer=5)
