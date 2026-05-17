"""Convert vessel skeleton edge-list parquet into nodes/edges graph tables."""

import pandas as pd

INPUT_EDGE_COLUMNS = [
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

NODE_COLUMNS = ["node_id", "channel", "vox_x", "vox_y", "vox_z", "x", "y", "z", "radius"]
EDGE_COLUMNS = [
    "edge_id",
    "channel",
    "src_node_id",
    "dst_node_id",
    "edge_length",
    "src_vox_x",
    "src_vox_y",
    "src_vox_z",
    "dst_vox_x",
    "dst_vox_y",
    "dst_vox_z",
]


def _empty_nodes():
    return pd.DataFrame(columns=NODE_COLUMNS)


def _empty_edges():
    return pd.DataFrame(columns=EDGE_COLUMNS)


def _validate_edge_table_columns(df):
    missing = [col for col in INPUT_EDGE_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(
            "Input vessel graph parquet is missing expected columns: "
            + ", ".join(missing)
        )


def build_nodes_table(edge_df):
    if edge_df.empty:
        return _empty_nodes()

    src_nodes = edge_df[
        ["channel", "src_vox_x", "src_vox_y", "src_vox_z", "src_x", "src_y", "src_z", "src_radius"]
    ].rename(
        columns={
            "src_vox_x": "vox_x",
            "src_vox_y": "vox_y",
            "src_vox_z": "vox_z",
            "src_x": "x",
            "src_y": "y",
            "src_z": "z",
            "src_radius": "radius",
        }
    )
    dst_nodes = edge_df[
        ["channel", "dst_vox_x", "dst_vox_y", "dst_vox_z", "dst_x", "dst_y", "dst_z", "dst_radius"]
    ].rename(
        columns={
            "dst_vox_x": "vox_x",
            "dst_vox_y": "vox_y",
            "dst_vox_z": "vox_z",
            "dst_x": "x",
            "dst_y": "y",
            "dst_z": "z",
            "dst_radius": "radius",
        }
    )

    nodes = pd.concat([src_nodes, dst_nodes], ignore_index=True)
    nodes = (
        nodes.groupby(["channel", "vox_x", "vox_y", "vox_z"], as_index=False)
        .agg({"x": "mean", "y": "mean", "z": "mean", "radius": "mean"})
        .sort_values(by=["channel", "vox_z", "vox_y", "vox_x"])
        .reset_index(drop=True)
    )
    nodes["node_id"] = nodes.index.astype("int64")
    return nodes[NODE_COLUMNS]


def build_edges_table(edge_df, nodes_df):
    if edge_df.empty:
        return _empty_edges()

    node_lookup = nodes_df[
        ["node_id", "channel", "vox_x", "vox_y", "vox_z"]
    ].rename(columns={"vox_x": "src_vox_x", "vox_y": "src_vox_y", "vox_z": "src_vox_z"})
    edges = edge_df[
        [
            "channel",
            "src_vox_x",
            "src_vox_y",
            "src_vox_z",
            "dst_vox_x",
            "dst_vox_y",
            "dst_vox_z",
            "edge_length",
        ]
    ].copy()

    edges = edges.merge(
        node_lookup.rename(columns={"node_id": "src_node_id"}),
        on=["channel", "src_vox_x", "src_vox_y", "src_vox_z"],
        how="left",
    )
    edges = edges.merge(
        node_lookup.rename(
            columns={
                "src_vox_x": "dst_vox_x",
                "src_vox_y": "dst_vox_y",
                "src_vox_z": "dst_vox_z",
                "src_node_id": "dst_node_id",
            }
        ),
        on=["channel", "dst_vox_x", "dst_vox_y", "dst_vox_z"],
        how="left",
    )

    if edges["src_node_id"].isna().any() or edges["dst_node_id"].isna().any():
        raise ValueError("Could not map all edge endpoints to node IDs.")

    edges["src_node_id"] = edges["src_node_id"].astype("int64")
    edges["dst_node_id"] = edges["dst_node_id"].astype("int64")
    edges = edges.sort_values(
        by=["channel", "src_node_id", "dst_node_id"]
    ).reset_index(drop=True)
    edges["edge_id"] = edges.index.astype("int64")
    return edges[EDGE_COLUMNS]


if __name__ == "__main__":
    edge_df = pd.read_parquet(snakemake.input.graph_parquet)
    _validate_edge_table_columns(edge_df)
    nodes_df = build_nodes_table(edge_df)
    edges_df = build_edges_table(edge_df, nodes_df)
    nodes_df.to_parquet(snakemake.output.nodes_parquet, index=False)
    edges_df.to_parquet(snakemake.output.edges_parquet, index=False)
