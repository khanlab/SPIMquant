from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd
import pytest


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


REPO_ROOT = _find_repo_root(Path(__file__).parent)
SCRIPT_PATH = (
    REPO_ROOT / "spimquant/workflow/scripts/convert_vessel_graph_to_nodes_edges.py"
)

spec = spec_from_file_location("convert_vessel_graph_to_nodes_edges", SCRIPT_PATH)
convert_mod = module_from_spec(spec)
spec.loader.exec_module(convert_mod)


def _sample_edges():
    return pd.DataFrame(
        [
            {
                "channel": 0,
                "src_vox_x": 1,
                "src_vox_y": 2,
                "src_vox_z": 3,
                "dst_vox_x": 2,
                "dst_vox_y": 2,
                "dst_vox_z": 3,
                "src_x": 1.0,
                "src_y": 2.0,
                "src_z": 3.0,
                "dst_x": 2.0,
                "dst_y": 2.0,
                "dst_z": 3.0,
                "src_radius": 0.5,
                "dst_radius": 0.6,
                "edge_length": 1.0,
            },
            {
                "channel": 0,
                "src_vox_x": 2,
                "src_vox_y": 2,
                "src_vox_z": 3,
                "dst_vox_x": 3,
                "dst_vox_y": 2,
                "dst_vox_z": 3,
                "src_x": 2.0,
                "src_y": 2.0,
                "src_z": 3.0,
                "dst_x": 3.0,
                "dst_y": 2.0,
                "dst_z": 3.0,
                "src_radius": 0.6,
                "dst_radius": 0.7,
                "edge_length": 1.0,
            },
        ]
    )


def test_build_nodes_and_edges_tables():
    edge_df = _sample_edges()
    convert_mod._validate_edge_table_columns(edge_df)

    nodes_df = convert_mod.build_nodes_table(edge_df)
    edges_df = convert_mod.build_edges_table(edge_df, nodes_df)

    assert list(nodes_df.columns) == convert_mod.NODE_COLUMNS
    assert list(edges_df.columns) == convert_mod.EDGE_COLUMNS
    assert len(nodes_df) == 3
    assert len(edges_df) == 2
    assert nodes_df["node_id"].is_unique
    assert edges_df["edge_id"].is_unique
    assert edges_df["src_node_id"].isin(nodes_df["node_id"]).all()
    assert edges_df["dst_node_id"].isin(nodes_df["node_id"]).all()


def test_empty_input_gives_empty_tables():
    edge_df = pd.DataFrame(columns=convert_mod.INPUT_EDGE_COLUMNS)
    nodes_df = convert_mod.build_nodes_table(edge_df)
    edges_df = convert_mod.build_edges_table(edge_df, nodes_df)

    assert list(nodes_df.columns) == convert_mod.NODE_COLUMNS
    assert list(edges_df.columns) == convert_mod.EDGE_COLUMNS
    assert nodes_df.empty
    assert edges_df.empty


def test_validate_columns_raises_for_missing():
    edge_df = _sample_edges().drop(columns=["src_x"])
    with pytest.raises(ValueError, match="missing expected columns"):
        convert_mod._validate_edge_table_columns(edge_df)


def test_streaming_nodes_edges_parquet_roundtrip(tmp_path):
    pytest.importorskip("pyarrow.parquet")
    edge_df = _sample_edges()
    graph_parquet = tmp_path / "graph.parquet"
    nodes_parquet = tmp_path / "nodes.parquet"
    edges_parquet = tmp_path / "edges.parquet"
    edge_df.to_parquet(graph_parquet, index=False)

    convert_mod.validate_input_parquet_columns(graph_parquet)
    nodes_df = convert_mod.build_nodes_table_from_parquet(graph_parquet, batch_size=1)
    nodes_df.to_parquet(nodes_parquet, index=False)
    convert_mod.write_edges_table_from_parquet(
        graph_parquet, nodes_df, edges_parquet, batch_size=1
    )

    out_nodes = pd.read_parquet(nodes_parquet)
    out_edges = pd.read_parquet(edges_parquet)
    assert list(out_nodes.columns) == convert_mod.NODE_COLUMNS
    assert list(out_edges.columns) == convert_mod.EDGE_COLUMNS
    assert len(out_nodes) == 3
    assert len(out_edges) == 2
