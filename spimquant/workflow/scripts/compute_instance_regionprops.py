"""Instance-level region properties: seeded watershed of a segmentation mask.

Alternative to compute_filtered_regionprops.py (connected components), for
masks where touching objects merge into large networks.  See instances.py for
the method.  Writes the same columns as the connected-components rule
(centroid and voxel count, named by ``regionprop_outputs``) so everything
downstream works unchanged, plus:

* ``seed_x/y/z`` -- position of the instance's seed (peak), physical units
* ``mean_intensity``, ``max_intensity`` -- of the corrected image
* ``peak_intensity`` -- smoothed corrected intensity at the seed
* ``elongation`` -- 1 for a sphere; large for neurite-like segments

and a JSON file with the parameters used, including how ``h`` was set.
"""

import json
from contextlib import ExitStack

import numpy as np

from instances import block_instances, block_noise_sd, h_from_noise, to_physical

# Fixed seed so the noise estimate (and hence h) is reproducible
_RANDOM_STATE = 42

# Gaussians are truncated at 4 sigma (scipy's default)
_TRUNCATE = 4.0


def _finite_or_none(v):
    """For JSON: round finite floats, map NaN/inf to null."""
    return round(float(v), 4) if np.isfinite(v) else None


def _block_grid(chunks):
    """Per-axis (start, size) of each chunk, from a dask ``.chunks`` tuple."""
    grid = []
    for axis_chunks in chunks:
        starts = np.concatenate([[0], np.cumsum(axis_chunks)[:-1]])
        grid.append(list(zip(starts.tolist(), list(axis_chunks))))
    return grid


def _estimate_h(params, img, mask, zooms_um, znimg_path, stat_level):
    """h from the noise in tissue-rich blocks, or the fallback."""
    import dask
    import dask.array as da
    from dogthresh import select_stat_blocks, tissue_fraction_per_block
    from zarrnii import ZarrNii

    fallback = float(params.h_fallback)
    # pick tissue-rich blocks from a coarse level, read in this process
    with dask.config.set(scheduler="threads"):
        level = stat_level
        while True:
            try:
                lo = ZarrNii.from_file(znimg_path, level=level)
                break
            except Exception:
                if level == 0:
                    raise
                level -= 1
        lowres = np.asarray(lo.data[0].compute())
    fractions = tissue_fraction_per_block(lowres, img.chunks, img.shape)
    del lowres
    blocks = select_stat_blocks(fractions, int(params.stat_blocks))
    if not blocks:
        return fallback, float("nan"), "fallback (no tissue blocks)", []

    depth = tuple(
        int(np.ceil(_TRUNCATE * (params.sigma_um + params.bg_sigma_um) / z))
        for z in zooms_um
    )
    oi = da.overlap.overlap(img, depth=depth, boundary=0)
    om = da.overlap.overlap(mask, depth=depth, boundary=0)
    sds = [
        dask.delayed(block_noise_sd)(
            oi.blocks[b],
            om.blocks[b],
            float(params.sigma_um),
            float(params.bg_sigma_um),
            zooms_um,
            int(params.stat_samples_per_block),
            _RANDOM_STATE + i,
        )
        for i, b in enumerate(blocks)
    ]
    sds = np.asarray(dask.compute(*sds), dtype=np.float64)
    h, sd, source = h_from_noise(
        sds, float(params.noise_k), fallback, float(params.noise_percentile)
    )
    return h, sd, source, sds.tolist()


def _block_task(img_block, mask_block, grid, depth, block_info=None, **kwargs):
    """``block_instances`` for one overlapped block, as a (1, 1, 1) object array."""
    idx = block_info[0]["chunk-location"]
    starts = [grid[a][i][0] for a, i in enumerate(idx)]
    sizes = [grid[a][i][1] for a, i in enumerate(idx)]
    out = np.empty((1, 1, 1), dtype=object)
    out[0, 0, 0] = block_instances(
        img_block,
        mask_block,
        core_lo=depth,
        core_hi=tuple(d + n for d, n in zip(depth, sizes)),
        offset=tuple(s - d for s, d in zip(starts, depth)),
        **kwargs,
    )
    return out


def _iter_results(client, tasks, batch):
    """Results of per-block Delayed ``tasks``, yielded as they become available.

    With a distributed client, all blocks are submitted at once and results
    are collected as they complete (the scheduler orders the work and frees
    each result once collected).  Otherwise blocks are computed in batches.
    """
    if client is not None:
        from dask.distributed import as_completed

        futures = client.compute(tasks)
        for future in as_completed(futures):
            res = future.result()
            future.release()
            yield res[0, 0, 0]
        return
    import dask

    for b0 in range(0, len(tasks), batch):
        for res in dask.compute(*tasks[b0 : b0 + batch]):
            yield res[0, 0, 0]


def main():
    import dask
    import dask.array as da
    import pyarrow as pa
    import pyarrow.parquet as pq
    from dask_setup import get_dask_client
    from zarrnii import ZarrNii

    params = snakemake.params

    with ExitStack() as stack:
        client = stack.enter_context(
            get_dask_client(
                snakemake.config["dask_scheduler"],
                snakemake.threads,
                threads_per_worker=int(params.dask_threads_per_worker),
                mem_mb=snakemake.resources.get("mem_mb"),
                reserve_mb=float(params.dask_reserve_mb),
            )
        )
        znimg = ZarrNii.from_file(snakemake.input.corrected)
        znmask = ZarrNii.from_file(snakemake.input.mask)
        if znimg.data.shape != znmask.data.shape:
            raise ValueError(
                f"Corrected image {znimg.data.shape} and mask {znmask.data.shape} "
                "differ in shape; both must be at segmentation_level"
            )
        zooms_um = [float(z) * 1000.0 for z in znimg.get_zooms()]  # mm -> um
        affine = np.asarray(znimg.affine.matrix)

        chunk = int(params.chunk_size)
        img = znimg.data[0].rechunk((chunk, chunk, chunk))
        mask = znmask.data[0].rechunk(img.chunks)

        # ------------------------------------------------------------ #
        # 1. Seed prominence h
        # ------------------------------------------------------------ #
        h_cfg = params.h
        if isinstance(h_cfg, str) and h_cfg.lower() == "auto":
            h, noise_sd, h_source, n_noise = _estimate_h(
                params,
                img,
                mask,
                zooms_um,
                snakemake.input.corrected,
                int(params.stat_level_offset),
            )
        else:
            h, noise_sd, h_source, n_noise = float(h_cfg), float("nan"), "config", []
        print(
            f"voxel size (um) {np.round(zooms_um, 3)}; h = {h:.3f} ({h_source}; "
            f"noise sd {noise_sd:.3f})",
            flush=True,
        )

        # ------------------------------------------------------------ #
        # 2. Instances, block by block; each owned by the block with its seed
        # ------------------------------------------------------------ #
        depth = tuple(int(np.ceil(float(params.depth_um) / z)) for z in zooms_um)
        oi = da.overlap.overlap(img, depth=depth, boundary=0)
        om = da.overlap.overlap(mask, depth=depth, boundary=0)
        # One graph for all blocks, split into per-block Delayed objects once.
        # (Building a delayed call per ``.blocks[idx]`` costs time proportional
        # to the whole graph, i.e. quadratic in the number of blocks: hours
        # for a whole brain.)
        results = da.map_blocks(
            _block_task,
            oi,
            om,
            dtype=object,
            chunks=tuple((1,) * n for n in img.numblocks),
            meta=np.empty((0, 0, 0), dtype=object),
            grid=_block_grid(img.chunks),
            depth=depth,
            zooms_um=zooms_um,
            sigma_um=float(params.sigma_um),
            min_dist_um=float(params.min_dist_um),
            h=h,
            min_volume_um3=float(params.min_volume_um3),
        )
        tasks = list(results.to_delayed().ravel())
        print(
            f"{len(tasks)} blocks of {img.chunksize}, overlap {depth} voxels",
            flush=True,
        )

        names = params.output_properties
        pos = names.get("centroid", "centroid")
        nvox = names.get("area", "area")
        columns = [f"{pos}_x", f"{pos}_y", f"{pos}_z", nvox]
        columns += ["seed_x", "seed_y", "seed_z"]
        columns += ["mean_intensity", "max_intensity", "peak_intensity", "elongation"]
        schema = pa.schema([pa.field(c, pa.float64()) for c in columns])

        n_instances = 0
        n_done = 0
        report_every = max(1, len(tasks) // 50)
        with pq.ParquetWriter(snakemake.output.regionprops_parquet, schema) as writer:
            for res in _iter_results(client, tasks, int(params.task_batch)):
                n_done += 1
                if res is not None:
                    cen = to_physical(res["centroid"], affine)
                    seed = to_physical(res["seed"], affine)
                    table = {
                        f"{pos}_x": cen[:, 0],
                        f"{pos}_y": cen[:, 1],
                        f"{pos}_z": cen[:, 2],
                        nvox: res["nvoxels"],
                        "seed_x": seed[:, 0],
                        "seed_y": seed[:, 1],
                        "seed_z": seed[:, 2],
                        "mean_intensity": res["mean_intensity"],
                        "max_intensity": res["max_intensity"],
                        "peak_intensity": res["peak_intensity"],
                        "elongation": res["elongation"],
                    }
                    writer.write_table(
                        pa.table(
                            {c: pa.array(table[c], type=pa.float64()) for c in columns}
                        )
                    )
                    n_instances += len(res["nvoxels"])
                if n_done % report_every == 0 or n_done == len(tasks):
                    print(
                        f"blocks {n_done}/{len(tasks)}: {n_instances:,} instances",
                        flush=True,
                    )

    stats = {
        "method": "seeded watershed of the mask from peaks of the smoothed image",
        "sigma_um": float(params.sigma_um),
        "min_dist_um": float(params.min_dist_um),
        "h": h,
        "h_source": h_source,
        "h_config": h_cfg,
        "noise_k": float(params.noise_k),
        "noise_sd": _finite_or_none(noise_sd),
        "noise_percentile": float(params.noise_percentile),
        "noise_sd_per_block": [_finite_or_none(v) for v in n_noise],
        "min_volume_um3": float(params.min_volume_um3),
        "depth_um": float(params.depth_um),
        "chunk_size": chunk,
        "n_instances": n_instances,
    }
    with open(snakemake.output.stats_json, "w") as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    main()
