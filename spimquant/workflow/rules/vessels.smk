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
        zarrnii_kwargs={"orientation": config["orientation"]},
        vesselfm_kwargs=lambda wildcards, input: {
            "chunk_size": (1, 128, 128, 128),
            "model_path": input.model_path,
        },
    output:
        mask=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="vesselfm",
                    suffix="mask.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            ),
            group_jobs=True,
        ),
    benchmark:
        bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="vesselfm",
            suffix="benchmark.tsv",
            **inputs["spim"].wildcards,
        )
    threads: 32
    resources:
        gpu=1,
        cpus_per_gpu=32,
        mem_mb=32000,
        runtime=lambda wildcards: max(1, int(200.0 / (3.0 ** float(wildcards.level)))),  # rough estimate, clamped to >=1
    script:
        "../scripts/vesselfm.py"


rule benchmark_run_vesselfm:
    """Benchmark VesselFM with different chunk sizes and thread counts.

    Runs VesselFM inference with varying chunk_size and threads to help
    identify optimal resource configurations. Results are written to a
    benchmark TSV file with timing and memory usage information.
    """
    input:
        spim=inputs["spim"].path,
        model_path="resources/models/vesselfm.pt",
    params:
        zarrnii_kwargs={"orientation": config["orientation"]},
        vesselfm_kwargs=lambda wildcards, input: {
            "chunk_size": (
                1,
                int(wildcards.chunk_size),
                int(wildcards.chunk_size),
                int(wildcards.chunk_size),
            ),
            "model_path": input.model_path,
        },
    output:
        mask=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="vesselfmBench",
                    res="{chunk_size}",
                    suffix="mask.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    benchmark:
        repeat(
            bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="vesselfmBench",
                res="{chunk_size}",
                nthreads="{thread_count}",
                suffix="benchmark.tsv",
                **inputs["spim"].wildcards,
            ),
            3,
        )
    threads: lambda wildcards: int(wildcards.thread_count)
    resources:
        gpu=1,
        cpus_per_gpu=lambda wildcards: int(wildcards.thread_count),
        mem_mb=32000,
        runtime=lambda wildcards: max(1, int(200.0 / (3.0 ** float(wildcards.level)))),
    script:
        "../scripts/vesselfm.py"


rule all_benchmark_vesselfm:
    """Aggregate benchmarks of VesselFM across chunk sizes and thread counts."""
    input:
        inputs["spim"].expand(
            bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="vesselfmBench",
                res="{chunk_size}",
                nthreads="{thread_count}",
                suffix="benchmark.tsv",
                **inputs["spim"].wildcards,
            ),
            stain=config["stains_for_vessels"],
            level=[config["segmentation_level"]],
            chunk_size=config["vesselfm_benchmark"]["chunk_sizes"],
            thread_count=config["vesselfm_benchmark"]["thread_counts"],
        ),


rule fieldfrac_vessels:
    """Calculate field fraction from binary mask.
    
    Computes the fraction of brain tissue occupied by the vessels.
    Note: This is a separate rule from `fieldfrac` to allow the 
    dags groups to be disjoint.
    
    """
    input:
        mask=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        hires_level=config["segmentation_level"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        fieldfrac_nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc,vesselfm}",
            suffix="fieldfrac.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/fieldfrac.py"


rule signed_distance_transform:
    """Compute signed distance transform from a binary mask.

    Applies the chamfer distance transform (distance_transform_cdt from scipy)
    to a binary mask using dask map_overlap for chunked, parallel processing.
    The output is a signed distance transform where positive values indicate
    the interior and negative values indicate the exterior of the mask.
    """
    input:
        mask=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="vesselfm",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        overlap_depth=32,
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        dist=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="{desc}",
                    suffix="dist.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            ),
            group_jobs=True,
        ),
    threads: 32
    resources:
        mem_mb=64000,
        disk_mb=2097152,
        runtime=180,
    script:
        "../scripts/signed_distance_transform.py"
