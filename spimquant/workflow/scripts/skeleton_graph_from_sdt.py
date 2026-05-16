"""Build a sparse vessel skeleton graph from skeleton mask + signed distance map.

This script reads a vessel skeleton OME-Zarr mask and a signed distance transform
(SDT) OME-Zarr map, then builds a sparse graph representation of the skeleton
using skan in chunked dask blocks. Node coordinates are converted to physical
space using OME-Zarr scale/translation metadata. Vessel radius at each skeleton
node is derived from SDT values sampled at skeleton voxels.
"""

import warnings

import numpy as np
import pandas as pd
from dask import compute, delayed
from dask.diagnostics import ProgressBar
from skan import Skeleton

EDGE_COLUMNS = [
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

AXIS_C = 0
AXIS_Z = 1
AXIS_Y = 2
AXIS_X = 3
EXPECTED_DIMS_CZYX = ("c", "z", "y", "x")
EXPECTED_DIMS_CXYZ = ("c", "x", "y", "z")


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


def _as_czyx_darr(zn_arr, input_name):
    """Return zarrnii `.darr` as (c, z, y, x) from either czyx or cxyz `.dims`.

    Expects an object exposing `dims` metadata and a dask-backed `darr` array.
    When dims are (c, x, y, z), transposes with (0, 3, 2, 1) to normalize to
    (c, z, y, x) for chunk processing.
    """
    dims = tuple(zn_arr.dims)
    if dims == EXPECTED_DIMS_CZYX:
        return zn_arr.darr
    if dims == EXPECTED_DIMS_CXYZ:
        # zarrnii can expose NIfTI-like arrays as (c, x, y, z); transpose to czyx.
        return zn_arr.darr.transpose(0, 3, 2, 1)
    raise ValueError(
        f"Expected dims {EXPECTED_DIMS_CZYX} or {EXPECTED_DIMS_CXYZ} for "
        f"{input_name}, got {dims}."
    )


def _has_neighbor_pair_26(binary_zyx):
    """Return True if any foreground voxel has a 26-connected neighbor.

    This pre-check filters degenerate skeleton chunks made only of isolated
    points, which do not contribute graph edges and can make skan fail.
    """
    shape = binary_zyx.shape
    for dz in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dz == 0 and dy == 0 and dx == 0:
                    continue
                src = binary_zyx[
                    max(0, -dz) : shape[0] - max(0, dz),
                    max(0, -dy) : shape[1] - max(0, dy),
                    max(0, -dx) : shape[2] - max(0, dx),
                ]
                dst = binary_zyx[
                    max(0, dz) : shape[0] - max(0, -dz),
                    max(0, dy) : shape[1] - max(0, -dy),
                    max(0, dx) : shape[2] - max(0, -dx),
                ]
                if np.any(src & dst):
                    return True
    return False


def _aggregate_block_tables(block_tables):
    """Aggregate per-chunk edge tables into one deduplicated DataFrame."""
    non_empty = [
        df for df in block_tables if isinstance(df, pd.DataFrame) and not df.empty
    ]

    if not non_empty:
        return pd.DataFrame(columns=EDGE_COLUMNS)

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
    return out_df


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

    # Channel axis has no overlap; spatial axes can have overlap.
    core_start_local = [0, 0, 0, 0]
    for d in (AXIS_Z, AXIS_Y, AXIS_X):
        if chunk_starts[d] > 0:
            core_start_local[d] = overlap_depth

    # Global origin of this overlapped block (before the chunk's core starts).
    # For interior chunks, map_overlap prepends `overlap_depth` voxels.
    overlap_block_starts = [chunk_starts[d] - core_start_local[d] for d in range(4)]

    # End indices include channel axis (with no overlap offset by design).
    core_end_local = [core_start_local[d] + chunk_sizes[d] for d in range(4)]

    rows = []

    for c in range(skel_block.shape[0]):
        skel_c = skel_block[c] > 0
        if not np.any(skel_c):
            continue

        skel_u8 = skel_c.astype(np.uint8)
        if not _has_neighbor_pair_26(skel_u8):
            # Degenerate chunks with isolated points have no valid graph edges.
            continue
        skeleton = Skeleton(skel_u8)
        coords_zyx = np.asarray(skeleton.coordinates, dtype=np.float64)

        if coords_zyx.shape[0] == 0:
            continue

        # In-core mask to avoid chunk-duplicate node ownership.
        in_core = (
            (coords_zyx[:, 0] >= core_start_local[AXIS_Z])
            & (coords_zyx[:, 0] < core_end_local[AXIS_Z])
            & (coords_zyx[:, 1] >= core_start_local[AXIS_Y])
            & (coords_zyx[:, 1] < core_end_local[AXIS_Y])
            & (coords_zyx[:, 2] >= core_start_local[AXIS_X])
            & (coords_zyx[:, 2] < core_end_local[AXIS_X])
        )

        # Convert local chunk coords -> global voxel coords.
        global_zyx = np.empty_like(coords_zyx)
        global_zyx[:, 0] = overlap_block_starts[AXIS_Z] + coords_zyx[:, 0]
        global_zyx[:, 1] = overlap_block_starts[AXIS_Y] + coords_zyx[:, 1]
        global_zyx[:, 2] = overlap_block_starts[AXIS_X] + coords_zyx[:, 2]

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
            # Convert to native Python-int tuples so edge-key equality/hashing is
            # consistent across pandas/parquet serialization boundaries.
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
    from dask_setup import get_dask_client
    from zarrnii import ZarrNii

    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):
        # Overlap keeps local skeleton connectivity across chunk boundaries.
        # vessels.smk sets overlap_depth=32 for SDT, skeletonization, and this rule.
        overlap_depth = int(snakemake.params.overlap_depth)

        zn_skel = ZarrNii.from_file(snakemake.input.skeleton)
        zn_sdt = ZarrNii.from_file(snakemake.input.sdt)

        skel_darr = _as_czyx_darr(zn_skel, "vessel skeleton graph extraction")
        sdt_darr = _as_czyx_darr(zn_sdt, "vessel SDT graph extraction")

        if skel_darr.shape != sdt_darr.shape:
            raise ValueError(
                "Skeleton and SDT arrays must have identical shape, got "
                f"{skel_darr.shape} vs {sdt_darr.shape}."
            )

        scale_xyz = np.array(
            [
                float(zn_skel.scale.get("x", 1.0)),
                float(zn_skel.scale.get("y", 1.0)),
                float(zn_skel.scale.get("z", 1.0)),
            ],
            dtype=np.float64,
        )
        translation = getattr(zn_skel, "translation", {})
        translation_xyz = np.array(
            [
                float(translation.get("x", 0.0)),
                float(translation.get("y", 0.0)),
                float(translation.get("z", 0.0)),
            ],
            dtype=np.float64,
        )

        # Keys are axis indices for dims (c, z, y, x), respectively.
        overlap_per_axis = {
            AXIS_C: 0,
            AXIS_Z: overlap_depth,
            AXIS_Y: overlap_depth,
            AXIS_X: overlap_depth,
        }
        skel_overlap = skel_darr.map_overlap(
            lambda x: x,
            depth=overlap_per_axis,
            boundary=0,
            trim=False,
            dtype=skel_darr.dtype,
        )
        sdt_overlap = sdt_darr.map_overlap(
            lambda x: x,
            depth=overlap_per_axis,
            boundary=0,
            trim=False,
            dtype=sdt_darr.dtype,
        )

        affine_obj = getattr(zn_skel, "affine", None)
        affine_matrix = None
        if affine_obj is not None:
            # zarrnii affine is usually an object with `.matrix`; if not, try array-like.
            affine_matrix = np.asarray(getattr(affine_obj, "matrix", affine_obj))
            if affine_matrix.shape != (4, 4):
                warnings.warn(
                    "Invalid affine shape in skeleton OME-Zarr metadata; falling back "
                    "to scale/translation-based coordinate conversion.",
                    stacklevel=2,
                )
                affine_matrix = None

        skel_blocks = skel_overlap.to_delayed()
        sdt_blocks = sdt_overlap.to_delayed()

        delayed_tables = []
        for chunk_idx in np.ndindex(*skel_overlap.numblocks):
            skel_block = skel_blocks[chunk_idx]
            sdt_block = sdt_blocks[chunk_idx]

            delayed_tables.append(
                delayed(_process_chunk)(
                    skel_block,
                    sdt_block,
                    chunk_idx,
                    skel_darr.chunks,
                    overlap_depth,
                    affine_matrix,
                    scale_xyz,
                    translation_xyz,
                )
            )

        with ProgressBar():
            block_tables = list(compute(*delayed_tables))

        out_df = _aggregate_block_tables(block_tables)
        out_df.to_parquet(snakemake.output.graph_parquet, index=False)
