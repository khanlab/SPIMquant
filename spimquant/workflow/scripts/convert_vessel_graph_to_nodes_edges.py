"""Convert vessel skeleton edge-list parquet into graph-package tables.

Expected Snakemake I/O
----------------------
Input:
    snakemake.input.graph_parquet
Outputs:
    snakemake.output.nodes_parquet
    snakemake.output.edges_parquet
"""

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

NODE_COLUMNS = [
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
PARQUET_BATCH_SIZE = 200_000


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


def _validate_input_parquet_columns(graph_parquet):
    """Validate required vessel graph columns from parquet schema."""
    import pyarrow.parquet as pq

    col_names = set(pq.ParquetFile(graph_parquet).schema.names)
    missing = [col for col in INPUT_EDGE_COLUMNS if col not in col_names]
    if missing:
        raise ValueError(
            "Input vessel graph parquet is missing expected columns: "
            + ", ".join(missing)
        )


def write_nodes_table_from_parquet(
    graph_parquet, nodes_parquet, batch_size=PARQUET_BATCH_SIZE
):
    """Build nodes table with Dask and write to parquet in bounded memory."""
    import dask.dataframe as dd
    import pyarrow as pa
    import pyarrow.parquet as pq

    source_ddf = dd.read_parquet(
        graph_parquet,
        columns=[
            "channel",
            "src_vox_x",
            "src_vox_y",
            "src_vox_z",
            "src_x",
            "src_y",
            "src_z",
            "src_radius",
            "dst_vox_x",
            "dst_vox_y",
            "dst_vox_z",
            "dst_x",
            "dst_y",
            "dst_z",
            "dst_radius",
        ],
        split_row_groups=True,
    )
    src_nodes = source_ddf[
        [
            "channel",
            "src_vox_x",
            "src_vox_y",
            "src_vox_z",
            "src_x",
            "src_y",
            "src_z",
            "src_radius",
        ]
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
    dst_nodes = source_ddf[
        [
            "channel",
            "dst_vox_x",
            "dst_vox_y",
            "dst_vox_z",
            "dst_x",
            "dst_y",
            "dst_z",
            "dst_radius",
        ]
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
    nodes_ddf = dd.concat([src_nodes, dst_nodes], interleave_partitions=True)
    nodes_ddf = (
        nodes_ddf.groupby(["channel", "vox_x", "vox_y", "vox_z"])[
            ["x", "y", "z", "radius"]
        ]
        .mean()
        .reset_index()
    )

    node_id_offset = 0
    writer = None
    try:
        for delayed_part in nodes_ddf.to_delayed():
            nodes_df = delayed_part.compute()
            if nodes_df.empty:
                continue
            nodes_df = nodes_df.sort_values(
                by=["channel", "vox_z", "vox_y", "vox_x"]
            ).reset_index(drop=True)
            n_rows = len(nodes_df)
            nodes_df["node_id"] = pd.RangeIndex(
                node_id_offset, node_id_offset + n_rows
            ).astype("int64")
            node_id_offset += n_rows
            node_table = pa.Table.from_pandas(
                nodes_df[NODE_COLUMNS], preserve_index=False
            )
            if writer is None:
                writer = pq.ParquetWriter(nodes_parquet, node_table.schema)
            writer.write_table(node_table)
    finally:
        if writer is not None:
            writer.close()

    if node_id_offset == 0:
        _empty_nodes().to_parquet(nodes_parquet, index=False)


def build_nodes_table(edge_df):
    if edge_df.empty:
        return _empty_nodes()

    src_nodes = edge_df[
        [
            "channel",
            "src_vox_x",
            "src_vox_y",
            "src_vox_z",
            "src_x",
            "src_y",
            "src_z",
            "src_radius",
        ]
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
        [
            "channel",
            "dst_vox_x",
            "dst_vox_y",
            "dst_vox_z",
            "dst_x",
            "dst_y",
            "dst_z",
            "dst_radius",
        ]
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

    node_lookup = nodes_df[["node_id", "channel", "vox_x", "vox_y", "vox_z"]].rename(
        columns={"vox_x": "src_vox_x", "vox_y": "src_vox_y", "vox_z": "src_vox_z"}
    )
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
                "node_id": "dst_node_id",
            }
        ),
        on=["channel", "dst_vox_x", "dst_vox_y", "dst_vox_z"],
        how="left",
    )

    if edges["src_node_id"].isna().any() or edges["dst_node_id"].isna().any():
        raise ValueError("Could not map all edge endpoints to node IDs.")

    edges["src_node_id"] = edges["src_node_id"].astype("int64")
    edges["dst_node_id"] = edges["dst_node_id"].astype("int64")
    edges = edges.sort_values(by=["channel", "src_node_id", "dst_node_id"]).reset_index(
        drop=True
    )
    edges["edge_id"] = edges.index.astype("int64")
    return edges[EDGE_COLUMNS]


def write_edges_table_from_parquet(
    graph_parquet, nodes_parquet, edges_parquet, batch_size=PARQUET_BATCH_SIZE
):
    """Map edges to node IDs with Dask and write output parquet."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    import dask.dataframe as dd

    edges_ddf = dd.read_parquet(
        graph_parquet,
        columns=[
            "channel",
            "src_vox_x",
            "src_vox_y",
            "src_vox_z",
            "dst_vox_x",
            "dst_vox_y",
            "dst_vox_z",
            "edge_length",
        ],
        split_row_groups=True,
    )
    nodes_lookup = dd.read_parquet(
        nodes_parquet,
        columns=["node_id", "channel", "vox_x", "vox_y", "vox_z"],
        split_row_groups=True,
    )
    src_lookup = nodes_lookup.rename(
        columns={
            "vox_x": "src_vox_x",
            "vox_y": "src_vox_y",
            "vox_z": "src_vox_z",
            "node_id": "src_node_id",
        }
    )
    dst_lookup = nodes_lookup.rename(
        columns={
            "vox_x": "dst_vox_x",
            "vox_y": "dst_vox_y",
            "vox_z": "dst_vox_z",
            "node_id": "dst_node_id",
        }
    )

    edges_with_nodes = edges_ddf.merge(
        src_lookup,
        on=["channel", "src_vox_x", "src_vox_y", "src_vox_z"],
        how="left",
    ).merge(
        dst_lookup,
        on=["channel", "dst_vox_x", "dst_vox_y", "dst_vox_z"],
        how="left",
    )

    import dask

    missing_src, missing_dst = dask.compute(
        edges_with_nodes["src_node_id"].isna().any(),
        edges_with_nodes["dst_node_id"].isna().any(),
    )
    if missing_src or missing_dst:
        raise ValueError("Could not map all edge endpoints to node IDs.")

    edge_id_offset = 0
    writer = None
    try:
        for delayed_part in edges_with_nodes.to_delayed():
            edge_df = delayed_part.compute()
            if edge_df.empty:
                continue
            edge_df = edge_df.reset_index(drop=True)
            n_rows = len(edge_df)
            edge_df["edge_id"] = pd.RangeIndex(
                edge_id_offset, edge_id_offset + n_rows
            ).astype("int64")
            edge_id_offset += n_rows
            edge_df["src_node_id"] = edge_df["src_node_id"].astype("int64")
            edge_df["dst_node_id"] = edge_df["dst_node_id"].astype("int64")
            edge_table = pa.Table.from_pandas(
                edge_df[EDGE_COLUMNS], preserve_index=False
            )
            if writer is None:
                writer = pq.ParquetWriter(edges_parquet, edge_table.schema)
            writer.write_table(edge_table)
    finally:
        if writer is not None:
            writer.close()

    if edge_id_offset == 0:
        _empty_edges().to_parquet(edges_parquet, index=False)


if __name__ == "__main__":
    from dask_setup import get_dask_client

    graph_parquet = snakemake.input.graph_parquet
    _validate_input_parquet_columns(graph_parquet)
    scheduler = snakemake.config.get("dask_scheduler", "threads")
    with get_dask_client(scheduler, snakemake.threads):
        write_nodes_table_from_parquet(
            graph_parquet,
            snakemake.output.nodes_parquet,
        )
        write_edges_table_from_parquet(
            graph_parquet,
            snakemake.output.nodes_parquet,
            snakemake.output.edges_parquet,
        )
