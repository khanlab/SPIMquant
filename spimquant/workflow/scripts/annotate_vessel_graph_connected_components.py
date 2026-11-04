"""Annotate vessel graph nodes with ranked connected-component labels.

Expected Snakemake I/O
----------------------
Inputs:
    snakemake.input.nodes_parquet
    snakemake.input.edges_parquet
Output:
    snakemake.output.nodes_parquet
"""

import pandas as pd

COMPONENT_COLUMN = "component_label"


def annotate_nodes_with_connected_components(nodes_df, edges_df):
    """Annotate each node with a ranked connected-component label.

    Components are ranked by descending size (largest first). Ties are resolved
    deterministically by the minimum node_id in each component.
    Labels start at 1.
    """
    annotated = nodes_df.copy()

    if annotated.empty:
        annotated[COMPONENT_COLUMN] = pd.Series(dtype="int64")
        return annotated

    if "node_id" not in annotated.columns:
        raise ValueError("nodes table must contain a node_id column")

    node_ids = annotated["node_id"].astype("int64").tolist()
    node_to_idx = {node_id: idx for idx, node_id in enumerate(node_ids)}

    parent = list(range(len(node_ids)))
    rank = [0] * len(node_ids)

    def find(i):
        """Return root index with path compression for amortized near-constant lookups."""
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        """Merge two sets using union-by-rank to keep trees shallow."""
        root_i = find(i)
        root_j = find(j)
        if root_i == root_j:
            return
        if rank[root_i] < rank[root_j]:
            parent[root_i] = root_j
        elif rank[root_i] > rank[root_j]:
            parent[root_j] = root_i
        else:
            parent[root_j] = root_i
            rank[root_i] += 1

    if not edges_df.empty:
        required_edge_cols = {"src_node_id", "dst_node_id"}
        missing = required_edge_cols.difference(edges_df.columns)
        if missing:
            raise ValueError(
                "edges table is missing expected columns: " + ", ".join(sorted(missing))
            )

        src_nodes = edges_df["src_node_id"].astype("int64")
        dst_nodes = edges_df["dst_node_id"].astype("int64")

        missing_node_ids = sorted(
            (set(src_nodes.unique()) | set(dst_nodes.unique())).difference(node_to_idx)
        )
        if missing_node_ids:
            raise ValueError(
                "edges reference node_ids not present in nodes table: "
                + ", ".join(map(str, missing_node_ids[:10]))
            )

        for src_id, dst_id in zip(src_nodes.tolist(), dst_nodes.tolist()):
            union(node_to_idx[src_id], node_to_idx[dst_id])

    components_by_root = {}
    for idx, node_id in enumerate(node_ids):
        root = find(idx)
        components_by_root.setdefault(root, []).append(node_id)

    ranked_components = sorted(
        components_by_root.values(),
        key=lambda comp_node_ids: (-len(comp_node_ids), min(comp_node_ids)),
    )

    component_label_by_node_id = {}
    for label, component_node_ids in enumerate(ranked_components, start=1):
        for node_id in component_node_ids:
            component_label_by_node_id[node_id] = label

    annotated[COMPONENT_COLUMN] = (
        annotated["node_id"]
        .astype("int64")
        .map(component_label_by_node_id)
        .astype("int64")
    )
    return annotated


if __name__ == "__main__":
    nodes_df = pd.read_parquet(snakemake.input.nodes_parquet)
    edges_df = pd.read_parquet(snakemake.input.edges_parquet)
    annotated_nodes = annotate_nodes_with_connected_components(nodes_df, edges_df)
    annotated_nodes.to_parquet(snakemake.output.nodes_parquet, index=False)
