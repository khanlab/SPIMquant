rule import_vesselfm_model:
    input:
        model=storage(config["models"]["vesselfm"]),
    output:
        "resources/models/vesselfm.pt",
    localrule: True
    shell:
        "cp {input} {output}"


rule run_vesselfm:
    """Segment vessels with vesselFM, via the shared tiled-inference engine.

    Overlapping 128^3 tiles (--vessel_stride) are merged across block
    boundaries before the sigmoid probability is thresholded; tiles and blocks
    outside the low-res brain mask are skipped (see vesselfm.py).
    """
    input:
        spim=spim_input,
        spim_json=spim_json_input,
        model_path="resources/models/vesselfm.pt",
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["correction_level"],
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
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
    threads: 8 * config["vessel_n_gpus"] if gpu_available else 8
    resources:
        # vesselfm runs on the CPU when torch has no CUDA; request a GPU only
        # where it can be used (see --gpu)
        gpu=config["vessel_n_gpus"] if gpu_available else 0,
        # all GPUs on one node: --gpus=N alone lets SLURM spread them over nodes,
        # and the job's single task then cannot use them ("Invalid generic
        # resource (gres) specification")
        nodes=1,
        # one task that sees every GPU (see run_lantern_plaques)
        tasks_per_gpu=0,
        cpus_per_gpu=8 if gpu_available else 0,
        # Up to `threads` blocks in flight (~3 GB each at chunk 512 + halo,
        # level 0) plus inference_read_procs readers (~1.6 GB each): measured
        # 21.5 GB with 4 readers (AS40F2, 4x, level 0).
        mem_mb=scaled_mem_mb(24000 + 1700 * config["inference_read_procs"], 20),
        # Inference splits across GPUs, the pyramid write does not. Measured on
        # AS40F2 (4x, level 0, 183 GB/channel): 67 min inference on 2 GPUs
        # (~17 s/block on 1), 27 min write; ~1.5x headroom here.
        runtime=scaled_runtime(30, 0.15 + 1.0 / config["vessel_n_gpus"]),
    params:
        zarrnii_kwargs=zarrnii_in_kwargs,
        tile=config["vessel_tile"],
        stride=config["vessel_stride"],
        batch_size=config["vessel_batch_size"],
        chunk=config["vessel_chunk"],
        merge_mode=config["vessel_merge_mode"],
        threshold=config["vessel_threshold"],
        n_gpus=config["vessel_n_gpus"],
        read_procs=config["inference_read_procs"],
    script:
        "../scripts/vesselfm.py"


rule signed_distance_transform:
    """Compute signed distance transform from a binary mask (vessels or seg).

    Applies the chamfer distance transform (distance_transform_cdt from scipy)
    to a binary mask using dask map_overlap for chunked, parallel processing.
    The output is a signed distance transform computed as dt_outside - dt_inside,
    where negative values indicate the interior and positive values indicate
    the exterior of the mask.
    """
    input:
        mask=bids_oz_in(
            root=root,
            datatype="{datatype}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    output:
        # Full resolution and float32, so about as large as the input; only
        # sampled at object centroids (and skeletons), so temporary unless
        # keep_sdt is set.
        dist=(
            bids_oz_out(
                root=root,
                datatype="{datatype,seg|vessels}",
                stain="{stain}",
                level="{level}",
                desc="{desc}",
                suffix="dist.{ext}",
                **inputs["spim"].wildcards,
            )
            if config["keep_sdt"]
            else temp(
                bids_oz_out(
                    root=root,
                    datatype="{datatype,seg|vessels}",
                    stain="{stain}",
                    level="{level}",
                    desc="{desc}",
                    suffix="dist.{ext}",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=scaled_mem_mb(200000, 75),
        disk_mb=2097152,
        runtime=scaled_runtime(30, 0.5),
    params:
        overlap_depth=32,
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
        mem_mb=scaled_mem_mb(32000, 60),
        disk_mb=2097152,
        runtime=scaled_runtime(15, 0.1),
    params:
        overlap_depth=32,
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
    params:
        overlap_depth=32,
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
        nodes_parquet=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level="{level}",
            desc="{desc}+skeleton",
            suffix="nodes.parquet",
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
