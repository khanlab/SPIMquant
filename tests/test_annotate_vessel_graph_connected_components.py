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
    REPO_ROOT
    / "spimquant/workflow/scripts/annotate_vessel_graph_connected_components.py"
)

spec = spec_from_file_location(
    "annotate_vessel_graph_connected_components", SCRIPT_PATH
)
mod = module_from_spec(spec)
spec.loader.exec_module(mod)


def _sample_nodes(node_ids):
    return pd.DataFrame(
        {
            "node_id": node_ids,
            "channel": [0] * len(node_ids),
            "vox_x": node_ids,
            "vox_y": [0] * len(node_ids),
            "vox_z": [0] * len(node_ids),
            "x": [float(n) for n in node_ids],
            "y": [0.0] * len(node_ids),
            "z": [0.0] * len(node_ids),
            "radius": [1.0] * len(node_ids),
        }
    )


def test_connected_components_ranked_by_size():
    nodes_df = _sample_nodes([0, 1, 2, 3, 4, 5])
    edges_df = pd.DataFrame(
        {
            "edge_id": [0, 1, 2],
            "channel": [0, 0, 0],
            "src_node_id": [0, 1, 3],
            "dst_node_id": [1, 2, 4],
        }
    )

    out = mod.annotate_nodes_with_connected_components(nodes_df, edges_df)

    labels = out.set_index("node_id")[mod.COMPONENT_COLUMN].to_dict()
    assert labels[0] == labels[1] == labels[2] == 1
    assert labels[3] == labels[4] == 2
    assert labels[5] == 3


def test_component_tie_breaks_by_min_node_id():
    nodes_df = _sample_nodes([0, 1, 2, 3])
    edges_df = pd.DataFrame(
        {
            "edge_id": [0, 1],
            "channel": [0, 0],
            "src_node_id": [0, 2],
            "dst_node_id": [1, 3],
        }
    )

    out = mod.annotate_nodes_with_connected_components(nodes_df, edges_df)

    labels = out.set_index("node_id")[mod.COMPONENT_COLUMN].to_dict()
    assert labels[0] == labels[1] == 1
    assert labels[2] == labels[3] == 2


def test_empty_nodes_keeps_schema_and_is_empty():
    nodes_df = pd.DataFrame(
        columns=[
            "node_id",
            "channel",
            "vox_x",
            "vox_y",
            "vox_z",
            "x",
            "y",
            "z",
            "radius",
        ]
    )
    edges_df = pd.DataFrame(columns=["src_node_id", "dst_node_id"])

    out = mod.annotate_nodes_with_connected_components(nodes_df, edges_df)

    assert out.empty
    assert mod.COMPONENT_COLUMN in out.columns


def test_edges_with_unknown_node_raise_error():
    nodes_df = _sample_nodes([0, 1])
    edges_df = pd.DataFrame(
        {
            "src_node_id": [0],
            "dst_node_id": [99],
        }
    )

    with pytest.raises(ValueError, match="edges reference node_ids"):
        mod.annotate_nodes_with_connected_components(nodes_df, edges_df)
