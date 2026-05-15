"""Build a sparse vessel skeleton graph from skeleton mask + signed distance map.

This script reads a vessel skeleton OME-Zarr mask and a signed distance transform
(SDT) OME-Zarr map, then builds a sparse graph representation of the skeleton
using skan in chunked dask blocks. Node coordinates are converted to physical
space using OME-Zarr scale/translation metadata. Vessel radius at each skeleton
node is derived from SDT values sampled at skeleton voxels.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dask import compute, delayed
from dask.diagnostics import ProgressBar
from skan import Skeleton
from zarrnii import ZarrNii

from dask_setup import get_dask_client


def _coord_from_scale_translation(voxel_xyz, scale_xyz, translation_xyz):
    """Convert voxel coordinates to physical coordinates."""
    return (voxel_xyz * scale_xyz) + translation_xyz


def _coord_from_affine(voxel_xyz, affine_matrix):
    """Convert voxel coordinates (x, y, z) to physical coordinates using affine."""
    voxel_h = np.concatenate(
        [voxel_xyz, np.ones((voxel_xyz.shape[0], 1), dtype=np.float64)],
        axis=1,
    )
    return (voxel_h @ affine_matrix.T)[:, :3]


def _process_chunk(
    skel_block,
    sdt_block,
    chunk_index,
    chunks,
    overlap_depth,
    affine_matrix,
    scale_xyz,
    translation_xyz,
):
    """Extract sparse graph edges from one overlapped chunk."""
    if skel_block.size == 0:
        return pd.DataFrame()

    # Compute global core bounds for this chunk (C, Z, Y, X)
    chunk_starts = [sum(chunks[d][: chunk_index[d]]) for d in range(4)]
    chunk_sizes = [chunks[d][chunk_index[d]] for d in range(4)]

    core_start_local = [0, 0, 0, 0]
    for d in range(1, 4):
        if chunk_starts[d] > 0:
            core_start_local[d] = overlap_depth

    core_end_local = [
        core_start_local[d] + chunk_sizes[d]
        for d in range(4)
    ]

    rows = []

    for c in range(skel_block.shape[0]):
        skel_c = skel_block[c] > 0
        if not np.any(skel_c):
            continue

        skeleton = Skeleton(skel_c.astype(np.uint8))
        coords_zyx = np.asarray(skeleton.coordinates, dtype=np.float64)

        if coords_zyx.shape[0] == 0:
            continue

        # In-core mask to avoid chunk-duplicate node ownership.
        in_core = (
            (coords_zyx[:, 0] >= core_start_local[1])
            & (coords_zyx[:, 0] < core_end_local[1])
            & (coords_zyx[:, 1] >= core_start_local[2])
            & (coords_zyx[:, 1] < core_end_local[2])
            & (coords_zyx[:, 2] >= core_start_local[3])
            & (coords_zyx[:, 2] < core_end_local[3])
        )

        # Convert local chunk coords -> global voxel coords.
        global_zyx = np.empty_like(coords_zyx)
        global_zyx[:, 0] = chunk_starts[1] + (coords_zyx[:, 0] - core_start_local[1])
        global_zyx[:, 1] = chunk_starts[2] + (coords_zyx[:, 1] - core_start_local[2])
        global_zyx[:, 2] = chunk_starts[3] + (coords_zyx[:, 2] - core_start_local[3])

        global_xyz = global_zyx[:, [2, 1, 0]]
        if affine_matrix is not None:
            phys_xyz = _coord_from_affine(global_xyz, affine_matrix)
        else:
            phys_xyz = _coord_from_scale_translation(
                global_xyz, scale_xyz, translation_xyz
            )

        # Radius from signed distance transform at skeleton voxels.
        local_zyx_idx = np.rint(coords_zyx).astype(np.int64)
        radius_vals = np.abs(
            sdt_block[
                c,
                local_zyx_idx[:, 0],
                local_zyx_idx[:, 1],
                local_zyx_idx[:, 2],
            ]
        ).astype(np.float64)

        graph = skeleton.graph.tocoo()

        for src, dst in zip(graph.row, graph.col):
            if src == dst:
                continue

            # Keep only one ownership direction: source node must be in core.
            if not in_core[src]:
                continue

            src_xyz_vox = global_xyz[src]
            dst_xyz_vox = global_xyz[dst]

            src_xyz_phys = phys_xyz[src]
            dst_xyz_phys = phys_xyz[dst]

            # Canonicalize edge orientation by voxel coordinate tuple to deduplicate.
            src_key = tuple(np.rint(src_xyz_vox).astype(np.int64).tolist())
            dst_key = tuple(np.rint(dst_xyz_vox).astype(np.int64).tolist())

            if src_key <= dst_key:
                p0_vox, p1_vox = src_xyz_vox, dst_xyz_vox
                p0_phys, p1_phys = src_xyz_phys, dst_xyz_phys
                r0, r1 = radius_vals[src], radius_vals[dst]
            else:
                p0_vox, p1_vox = dst_xyz_vox, src_xyz_vox
                p0_phys, p1_phys = dst_xyz_phys, src_xyz_phys
                r0, r1 = radius_vals[dst], radius_vals[src]

            edge_len = float(np.linalg.norm(p1_phys - p0_phys))

            rows.append(
                {
                    "channel": int(c),
                    "src_vox_x": int(np.rint(p0_vox[0])),
                    "src_vox_y": int(np.rint(p0_vox[1])),
                    "src_vox_z": int(np.rint(p0_vox[2])),
                    "dst_vox_x": int(np.rint(p1_vox[0])),
                    "dst_vox_y": int(np.rint(p1_vox[1])),
                    "dst_vox_z": int(np.rint(p1_vox[2])),
                    "src_x": float(p0_phys[0]),
                    "src_y": float(p0_phys[1]),
                    "src_z": float(p0_phys[2]),
                    "dst_x": float(p1_phys[0]),
                    "dst_y": float(p1_phys[1]),
                    "dst_z": float(p1_phys[2]),
                    "src_radius": float(r0),
                    "dst_radius": float(r1),
                    "edge_length": edge_len,
                }
            )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    # Remove duplicates that can still arise from overlap ownership in rare cases.
    dedup_cols = [
        "channel",
        "src_vox_x",
        "src_vox_y",
        "src_vox_z",
        "dst_vox_x",
        "dst_vox_y",
        "dst_vox_z",
    ]
    return df.drop_duplicates(subset=dedup_cols, keep="first")


if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        overlap_depth = int(snakemake.params.overlap_depth)

        zn_skel = ZarrNii.from_file(snakemake.input.skeleton)
        zn_sdt = ZarrNii.from_file(snakemake.input.sdt)

        expected_dims = ("c", "z", "y", "x")
        if tuple(zn_skel.dims) != expected_dims:
            raise ValueError(
                f"Expected dims {expected_dims} for vessel skeleton graph extraction, "
                f"got {tuple(zn_skel.dims)}."
            )
        if tuple(zn_sdt.dims) != expected_dims:
            raise ValueError(
                f"Expected dims {expected_dims} for vessel SDT graph extraction, "
                f"got {tuple(zn_sdt.dims)}."
            )
        if zn_skel.darr.shape != zn_sdt.darr.shape:
            raise ValueError(
                "Skeleton and SDT arrays must have identical shape, got "
                f"{zn_skel.darr.shape} vs {zn_sdt.darr.shape}."
            )

        scale_xyz = np.array(
            [
                float(zn_skel.scale.get("x", 1.0)),
                float(zn_skel.scale.get("y", 1.0)),
                float(zn_skel.scale.get("z", 1.0)),
            ],
            dtype=np.float64,
        )
        translation = getattr(zn_skel, "translation", {}) or {}
        translation_xyz = np.array(
            [
                float(translation.get("x", 0.0)),
                float(translation.get("y", 0.0)),
                float(translation.get("z", 0.0)),
            ],
            dtype=np.float64,
        )

        depth = {0: 0, 1: overlap_depth, 2: overlap_depth, 3: overlap_depth}
        skel_overlap = zn_skel.darr.map_overlap(
            lambda x: x,
            depth=depth,
            boundary=0,
            trim=False,
            dtype=zn_skel.darr.dtype,
        )
        sdt_overlap = zn_sdt.darr.map_overlap(
            lambda x: x,
            depth=depth,
            boundary=0,
            trim=False,
            dtype=zn_sdt.darr.dtype,
        )

        affine_obj = getattr(zn_skel, "affine", None)
        affine_matrix = None
        if affine_obj is not None:
            affine_matrix = np.asarray(getattr(affine_obj, "matrix", affine_obj))
            if affine_matrix.shape != (4, 4):
                affine_matrix = None

        skel_blocks = skel_overlap.to_delayed()
        sdt_blocks = sdt_overlap.to_delayed()

        delayed_tables = []
        for block_idx in np.ndindex(*skel_overlap.numblocks):
            skel_block = skel_blocks[block_idx]
            sdt_block = sdt_blocks[block_idx]

            delayed_tables.append(
                delayed(_process_chunk)(
                    skel_block,
                    sdt_block,
                    block_idx,
                    zn_skel.darr.chunks,
                    overlap_depth,
                    affine_matrix,
                    scale_xyz,
                    translation_xyz,
                )
            )

        with ProgressBar():
            block_tables = list(compute(*delayed_tables))

        non_empty = [
            df for df in block_tables if isinstance(df, pd.DataFrame) and not df.empty
        ]

        if non_empty:
            out_df = pd.concat(non_empty, ignore_index=True)
            out_df = out_df.drop_duplicates(
                subset=[
                    "channel",
                    "src_vox_x",
                    "src_vox_y",
                    "src_vox_z",
                    "dst_vox_x",
                    "dst_vox_y",
                    "dst_vox_z",
                ],
                keep="first",
            )
            out_df = out_df.sort_values(
                by=[
                    "channel",
                    "src_vox_z",
                    "src_vox_y",
                    "src_vox_x",
                    "dst_vox_z",
                    "dst_vox_y",
                    "dst_vox_x",
                ]
            ).reset_index(drop=True)
        else:
            out_df = pd.DataFrame(
                columns=[
                    "channel",
                    "src_vox_x",
                    "src_vox_y",
                    "src_vox_z",
                    "dst_vox_x",
                    "dst_vox_y",
                    "dst_vox_z",
                    "src_x",
                    "src_y",
                    "src_z",
                    "dst_x",
                    "dst_y",
                    "dst_z",
                    "src_radius",
                    "dst_radius",
                    "edge_length",
                ]
            )

        out_df.to_parquet(snakemake.output.graph_parquet, index=False)
