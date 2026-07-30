rule import_vesselfm_model:
    input:
        model=storage(config["models"]["vesselfm"]),
    output:
        "resources/models/vesselfm.pt",
    localrule: True
    shell:
        "cp {input} {output}"


rule run_vesselfm:
    input:
        spim=inputs["spim"].path,
        model_path="resources/models/vesselfm.pt",
    params:
        zarrnii_kwargs=zarrnii_in_kwargs,
        vesselfm_kwargs=lambda wildcards, input: {
            "chunk_size": (1, 128, 128, 128),
            "model_path": input.model_path,
        },
    output:
        mask=bids_oz_out(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="vesselfm",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    threads: 8
    resources:
        gpu=1,
        cpus_per_gpu=8,
        mem_mb=256000,
        runtime=lambda wildcards: max(1, int(400.0 / (3.0 ** float(wildcards.level)))),  # rough estimate, clamped to >=1
    script:
        "../scripts/vesselfm.py"


rule signed_distance_transform:
    """Compute signed distance transform from a binary mask.

    Applies the chamfer distance transform (distance_transform_cdt from scipy)
    to a binary mask using dask map_overlap for chunked, parallel processing.
    The output is a signed distance transform computed as dt_outside - dt_inside,
    where negative values indicate the interior and positive values indicate
    the exterior of the mask.
    """
    input:
        mask=bids_oz_in(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="vesselfm",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    params:
        overlap_depth=32,
    output:
        dist=bids_oz_out(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="dist.{ext}",
            **inputs["spim"].wildcards,
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=720,
    script:
        "../scripts/signed_distance_transform.py"


rule skeletonize_vessels_mask:
    """Skeletonize a vessel mask using chunked overlap-aware processing."""
    input:
        mask=bids_oz_in(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    params:
        overlap_depth=32,
    output:
        mask=bids_oz_out(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=720,
    script:
        "../scripts/skeletonize_vessels_mask.py"


rule vessel_skeleton_graph:
    """Create a sparse vessel graph parquet from skeleton mask and SDT."""
    input:
        skeleton=bids_oz_in(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
        sdt=bids_oz_in(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="dist.{ext}",
            **inputs["spim"].wildcards,
        ),
    params:
        overlap_depth=32,
    output:
        graph_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="graph.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=360,
    script:
        "../scripts/skeleton_graph_from_sdt.py"


rule vessel_graph_to_nodes_edges:
    """Convert vessel skeleton edge-list into graph-friendly node/edge tables."""
    input:
        graph_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="graph.parquet",
            **inputs["spim"].wildcards,
        ),
    output:
        nodes_raw_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="nodes_raw.parquet",
            **inputs["spim"].wildcards,
        ),
        edges_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="edges.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 128
    resources:
        mem_mb=256000,
        runtime=360,
    script:
        "../scripts/convert_vessel_graph_to_nodes_edges.py"


rule vessel_graph_connected_components:
    """Annotate vessel graph nodes with ranked connected-component labels."""
    input:
        nodes_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="nodes_raw.parquet",
            **inputs["spim"].wildcards,
        ),
        edges_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="edges.parquet",
            **inputs["spim"].wildcards,
        ),
    output:
        nodes_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="nodes.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 16
    resources:
        mem_mb=64000,
        runtime=360,
    script:
        "../scripts/annotate_vessel_graph_connected_components.py"
