from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
import pytest
from scipy.ndimage import distance_transform_edt

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "spimquant/workflow/scripts/skeleton_graph_from_sdt.py"


spec = spec_from_file_location("skeleton_graph_from_sdt", SCRIPT_PATH)
skeleton_graph_mod = module_from_spec(spec)
spec.loader.exec_module(skeleton_graph_mod)


def _load_example_or_synthetic():
    """Load provided example masks when present; otherwise create synthetic masks."""
    vessels_nii = REPO_ROOT / "tests/example_vessels_mask.nii.gz"
    skeleton_nii = REPO_ROOT / "tests/example_skeleton.nii.gz"

    if vessels_nii.exists() and skeleton_nii.exists():
        nib = pytest.importorskip("nibabel")
        vessels = nib.load(vessels_nii).get_fdata().astype(np.uint8)
        skeleton = nib.load(skeleton_nii).get_fdata().astype(np.uint8)
        return vessels, skeleton

    shape = (64, 64, 64)
    vessels = np.zeros(shape, dtype=np.uint8)
    skeleton = np.zeros(shape, dtype=np.uint8)

    # Simple centerline-like structure and a slightly thicker vessel mask.
    skeleton[30:34, 20:44, 32] = 100
    skeleton[32, 32, 18:46] = 100

    vessels[29:35, 19:45, 31:34] = 100
    vessels[31:34, 31:34, 17:47] = 100

    return vessels, skeleton


def _signed_distance_from_mask(mask):
    binary = mask > 0
    inside = distance_transform_edt(binary)
    outside = distance_transform_edt(~binary)
    return (outside - inside).astype(np.float32)


def _run_chunked_extraction(skeleton_zyx, sdt_zyx, overlap_depth=4, spatial_chunk=24):
    skel_4d = skeleton_zyx[np.newaxis, ...]
    sdt_4d = sdt_zyx[np.newaxis, ...]

    skel_da = da.from_array(
        skel_4d,
        chunks=(1, spatial_chunk, spatial_chunk, spatial_chunk),
    )
    sdt_da = da.from_array(
        sdt_4d,
        chunks=(1, spatial_chunk, spatial_chunk, spatial_chunk),
    )

    overlap = {
        skeleton_graph_mod.AXIS_C: 0,
        skeleton_graph_mod.AXIS_Z: overlap_depth,
        skeleton_graph_mod.AXIS_Y: overlap_depth,
        skeleton_graph_mod.AXIS_X: overlap_depth,
    }

    skel_overlap = skel_da.map_overlap(
        lambda x: x,
        depth=overlap,
        boundary=0,
        trim=False,
        dtype=skel_da.dtype,
    )
    sdt_overlap = sdt_da.map_overlap(
        lambda x: x,
        depth=overlap,
        boundary=0,
        trim=False,
        dtype=sdt_da.dtype,
    )

    block_tables = []
    for idx in np.ndindex(*skel_overlap.numblocks):
        skel_block = skel_overlap.blocks[idx].compute()
        sdt_block = sdt_overlap.blocks[idx].compute()
        block_tables.append(
            skeleton_graph_mod._process_chunk(
                skel_block,
                sdt_block,
                idx,
                skel_da.chunks,
                overlap_depth,
                None,
                np.array([1.0, 1.0, 1.0], dtype=np.float64),
                np.array([0.0, 0.0, 0.0], dtype=np.float64),
            )
        )

    return skeleton_graph_mod._aggregate_block_tables(block_tables)


def test_aggregate_block_tables_handles_empty_entries():
    out = skeleton_graph_mod._aggregate_block_tables([None, "not-a-df", pd.DataFrame()])
    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == skeleton_graph_mod.EDGE_COLUMNS
    assert out.empty


def test_process_chunk_and_chunked_aggregation_generate_graph_without_duplicates():
    vessels, skeleton = _load_example_or_synthetic()
    sdt = _signed_distance_from_mask(vessels)

    out = _run_chunked_extraction(skeleton, sdt, overlap_depth=4, spatial_chunk=24)

    assert isinstance(out, pd.DataFrame)
    assert list(out.columns) == skeleton_graph_mod.EDGE_COLUMNS
    assert len(out) > 0

    dedup_cols = [
        "channel",
        "src_vox_x",
        "src_vox_y",
        "src_vox_z",
        "dst_vox_x",
        "dst_vox_y",
        "dst_vox_z",
    ]
    assert not out.duplicated(subset=dedup_cols).any()
    assert (out["src_radius"] >= 0).all()
    assert (out["dst_radius"] >= 0).all()
    assert (out["edge_length"] > 0).all()


def test_affine_coordinate_conversion_is_applied():
    skeleton = np.zeros((20, 20, 20), dtype=np.uint8)
    skeleton[10, 5:9, 10] = 100
    sdt = np.ones_like(skeleton, dtype=np.float32)

    skel_block = skeleton[np.newaxis, ...]
    sdt_block = sdt[np.newaxis, ...]

    affine = np.array(
        [
            [2.0, 0.0, 0.0, 10.0],
            [0.0, 3.0, 0.0, 20.0],
            [0.0, 0.0, 4.0, 30.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )

    out = skeleton_graph_mod._process_chunk(
        skel_block,
        sdt_block,
        (0, 0, 0, 0),
        ((1,), (20,), (20,), (20,)),
        0,
        affine,
        np.array([1.0, 1.0, 1.0]),
        np.array([0.0, 0.0, 0.0]),
    )

    assert len(out) > 0
    # Physical coords should include affine translation and scaling (not voxel coords).
    assert (out["src_x"] >= 10.0).all()
    assert (out["src_y"] >= 20.0).all()
    assert (out["src_z"] >= 30.0).all()
