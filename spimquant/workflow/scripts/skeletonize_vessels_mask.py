"""Skeletonize a vessel mask with dask map_overlap and scikit-image."""

import dask.array as da
import numpy as np
from dask.diagnostics import ProgressBar
from skimage.morphology import skeletonize
from zarrnii import ZarrNii

from dask_setup import get_dask_client

MASK_TRUE_VALUE = np.uint8(100)


if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        overlap_depth = snakemake.params.overlap_depth

        znimg = ZarrNii.from_file(snakemake.input.mask)
        expected_dims = ("c", "z", "y", "x")
        if tuple(znimg.dims) != expected_dims:
            raise ValueError(
                f"Expected dims {expected_dims} for vessel mask skeletonization, "
                f"got {tuple(znimg.dims)}."
            )

        def skeletonize_block(block):
            """Skeletonize a chunk with shape (C, Z, Y, X) and return same shape."""
            result = np.zeros(block.shape, dtype=np.uint8)
            for c in range(block.shape[0]):
                binary = block[c] > 0
                if np.any(binary):
                    # Preserve the project's mask convention of foreground=100.
                    result[c] = skeletonize(binary).astype(np.uint8) * MASK_TRUE_VALUE
            return result

        # Dims are validated as (C, Z, Y, X): no overlap on C, overlap on Z/Y/X.
        depth = {0: 0, 1: overlap_depth, 2: overlap_depth, 3: overlap_depth}
        skel_darr = da.map_overlap(
            skeletonize_block,
            znimg.darr,
            depth=depth,
            boundary=0,
            dtype=np.uint8,
        )

        znimg.darr = skel_darr

        with ProgressBar():
            # Match pyramid generation used by other mask-producing scripts.
            znimg.to_ome_zarr(
                snakemake.output.mask,
                match_scale_factors_from=snakemake.input.mask,
                **snakemake.config["zarrnii_out_kwargs"],
            )
