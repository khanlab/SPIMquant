"""Compute signed distance transform from a binary mask using dask map_overlap.

The signed distance transform assigns:
- Positive values to interior voxels (inside the mask), proportional to distance
  to the nearest background voxel.
- Negative values to exterior voxels (outside the mask), proportional to distance
  to the nearest foreground voxel.

Uses scipy.ndimage.distance_transform_cdt via dask.array.map_overlap for chunked,
parallel processing with spatial overlap padding to reduce boundary artifacts.
"""

if __name__ == "__main__":
    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=snakemake.threads // 2,
        threads_per_worker=2,
        memory_limit="auto",
        dashboard_address=":8788",
    )
    client = Client(cluster)
    print(cluster.dashboard_link)

    try:
        import dask.array as da
        import numpy as np
        from scipy.ndimage import distance_transform_cdt
        from zarrnii import ZarrNii

        overlap_depth = snakemake.params.overlap_depth

        znimg = ZarrNii.from_ome_zarr(
            snakemake.input.mask, **snakemake.params.zarrnii_kwargs
        )

        def signed_dt_block(block):
            """Compute signed distance transform for a single block.

            Expects block shape (C, Z, Y, X).  For each channel the chamfer
            distance transform is computed twice:
              - dt_inside : distance from each foreground voxel to the nearest
                background voxel (positive contribution inside the mask).
              - dt_outside: distance from each background voxel to the nearest
                foreground voxel (positive contribution outside the mask).

            The signed distance transform is dt_inside - dt_outside, giving
            positive values inside the mask and negative values outside.
            """
            result = np.zeros(block.shape, dtype=np.float32)
            for c in range(block.shape[0]):
                binary = block[c] > 0
                dt_inside = distance_transform_cdt(binary).astype(np.float32)
                dt_outside = distance_transform_cdt(~binary).astype(np.float32)
                result[c] = dt_inside - dt_outside
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

    finally:
        client.close()
        cluster.close()
