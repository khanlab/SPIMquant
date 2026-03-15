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
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="vesselfm",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    resources:
        gpu=1,
        cpus_per_gpu=32,
        mem_mb=64000,
        runtime=lambda wildcards: max(1, int(200.0 / (3.0 ** float(wildcards.level)))),  # rough estimate, clamped to >=1
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
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="vesselfm",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    params:
        overlap_depth=32,
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        dist=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="dist.ozx",
            **inputs["spim"].wildcards,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=360,
    script:
        "../scripts/signed_distance_transform.py"
