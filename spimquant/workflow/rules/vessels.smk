rule run_vesselfm:
    input:
        spim=inputs["spim"].path,
        model=remote(
            "https://huggingface.co/bwittmann/vesselFM/resolve/main/vesselFM_base.pt"
        ),
    params:
        zarrnii_kwargs={"orientation": config["orientation"]},
        vesselfm_kwargs=lambda wildcards, input: {
            "chunk_size": (1, 128, 128, 128),
            "model": input.model,
        },
    output:
        mask=directory(
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
    threads: 32
    group:
        "subj"
    resources:
        mem_mb=32000,
        runtime=lambda wildcards: max(1, int(200.0 / (3.0 ** float(wildcards.level)))),  # rough estimate, clamped to >=1
    script:
        "../scripts/vesselfm.py"
