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


def validate_input_parquet_columns(graph_parquet):
    import pyarrow.parquet as pq

    col_names = set(pq.ParquetFile(graph_parquet).schema.names)
    missing = [col for col in INPUT_EDGE_COLUMNS if col not in col_names]
    if missing:
        raise ValueError(
            "Input vessel graph parquet is missing expected columns: "
            + ", ".join(missing)
        )


def _accumulate_endpoint_nodes(edge_df, endpoint, node_stats):
    vox_x = edge_df[f"{endpoint}_vox_x"].astype("int64").to_numpy()
    vox_y = edge_df[f"{endpoint}_vox_y"].astype("int64").to_numpy()
    vox_z = edge_df[f"{endpoint}_vox_z"].astype("int64").to_numpy()
    channels = edge_df["channel"].astype("int64").to_numpy()
    x = edge_df[f"{endpoint}_x"].astype("float64").to_numpy()
    y = edge_df[f"{endpoint}_y"].astype("float64").to_numpy()
    z = edge_df[f"{endpoint}_z"].astype("float64").to_numpy()
    radius = edge_df[f"{endpoint}_radius"].astype("float64").to_numpy()

    for c, vx, vy, vz, px, py, pz, pr in zip(
        channels, vox_x, vox_y, vox_z, x, y, z, radius
    ):
        key = (int(c), int(vx), int(vy), int(vz))
        if key not in node_stats:
            node_stats[key] = [float(px), float(py), float(pz), float(pr), 1]
        else:
            node_stats[key][0] += float(px)
            node_stats[key][1] += float(py)
            node_stats[key][2] += float(pz)
            node_stats[key][3] += float(pr)
            node_stats[key][4] += 1


def build_nodes_table_from_parquet(graph_parquet, batch_size=PARQUET_BATCH_SIZE):
    import pyarrow.parquet as pq

    node_stats = {}
    endpoint_columns = [
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
    ]

    parquet_file = pq.ParquetFile(graph_parquet)
    for batch in parquet_file.iter_batches(
        batch_size=batch_size, columns=endpoint_columns
    ):
        batch_df = batch.to_pandas()
        if batch_df.empty:
            continue
        _accumulate_endpoint_nodes(batch_df, "src", node_stats)
        _accumulate_endpoint_nodes(batch_df, "dst", node_stats)

    if not node_stats:
        return _empty_nodes()

    nodes = pd.DataFrame(
        (
            {
                "channel": channel,
                "vox_x": vox_x,
                "vox_y": vox_y,
                "vox_z": vox_z,
                "x": x_sum / count,
                "y": y_sum / count,
                "z": z_sum / count,
                "radius": radius_sum / count,
            }
            for (channel, vox_x, vox_y, vox_z), (
                x_sum,
                y_sum,
                z_sum,
                radius_sum,
                count,
            ) in node_stats.items()
        )
    )
    nodes = nodes.sort_values(by=["channel", "vox_z", "vox_y", "vox_x"]).reset_index(
        drop=True
    )
    nodes["node_id"] = nodes.index.astype("int64")
    return nodes[NODE_COLUMNS]


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
    edges = edges.sort_values(
        by=["channel", "src_node_id", "dst_node_id"]
    ).reset_index(drop=True)
    edges["edge_id"] = edges.index.astype("int64")
    return edges[EDGE_COLUMNS]


def write_edges_table_from_parquet(
    graph_parquet, nodes_df, edges_parquet, batch_size=PARQUET_BATCH_SIZE
):
    import pyarrow as pa
    import pyarrow.parquet as pq

    if nodes_df.empty:
        _empty_edges().to_parquet(edges_parquet, index=False)
        return

    node_lookup = {
        (int(row.channel), int(row.vox_x), int(row.vox_y), int(row.vox_z)): int(
            row.node_id
        )
        for row in nodes_df.itertuples(index=False)
    }

    parquet_file = pq.ParquetFile(graph_parquet)
    columns = [
        "channel",
        "src_vox_x",
        "src_vox_y",
        "src_vox_z",
        "dst_vox_x",
        "dst_vox_y",
        "dst_vox_z",
        "edge_length",
    ]
    edge_id_offset = 0
    writer = None

    try:
        for batch in parquet_file.iter_batches(batch_size=batch_size, columns=columns):
            batch_df = batch.to_pandas()
            if batch_df.empty:
                continue

            channels = batch_df["channel"].astype("int64").to_numpy()
            src_vox_x = batch_df["src_vox_x"].astype("int64").to_numpy()
            src_vox_y = batch_df["src_vox_y"].astype("int64").to_numpy()
            src_vox_z = batch_df["src_vox_z"].astype("int64").to_numpy()
            dst_vox_x = batch_df["dst_vox_x"].astype("int64").to_numpy()
            dst_vox_y = batch_df["dst_vox_y"].astype("int64").to_numpy()
            dst_vox_z = batch_df["dst_vox_z"].astype("int64").to_numpy()

            try:
                src_node_id = [
                    node_lookup[(c, sx, sy, sz)]
                    for c, sx, sy, sz in zip(channels, src_vox_x, src_vox_y, src_vox_z)
                ]
                dst_node_id = [
                    node_lookup[(c, dx, dy, dz)]
                    for c, dx, dy, dz in zip(channels, dst_vox_x, dst_vox_y, dst_vox_z)
                ]
            except KeyError as err:
                raise ValueError("Could not map all edge endpoints to node IDs.") from err

            n_rows = len(batch_df)
            edge_df = pd.DataFrame(
                {
                    "edge_id": range(edge_id_offset, edge_id_offset + n_rows),
                    "channel": channels,
                    "src_node_id": src_node_id,
                    "dst_node_id": dst_node_id,
                    "edge_length": batch_df["edge_length"].to_numpy(),
                    "src_vox_x": src_vox_x,
                    "src_vox_y": src_vox_y,
                    "src_vox_z": src_vox_z,
                    "dst_vox_x": dst_vox_x,
                    "dst_vox_y": dst_vox_y,
                    "dst_vox_z": dst_vox_z,
                }
            )
            edge_id_offset += n_rows

            edge_table = pa.Table.from_pandas(edge_df[EDGE_COLUMNS], preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(edges_parquet, edge_table.schema)
            writer.write_table(edge_table)
    finally:
        if writer is not None:
            writer.close()

    if writer is None:
        _empty_edges().to_parquet(edges_parquet, index=False)


if __name__ == "__main__":
    graph_parquet = snakemake.input.graph_parquet
    validate_input_parquet_columns(graph_parquet)
    nodes_df = build_nodes_table_from_parquet(graph_parquet)
    nodes_df.to_parquet(snakemake.output.nodes_parquet, index=False)
    write_edges_table_from_parquet(
        graph_parquet,
        nodes_df,
        snakemake.output.edges_parquet,
    )
