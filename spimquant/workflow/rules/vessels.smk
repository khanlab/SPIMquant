rule run_vesselfm:
    input:
        spim=inputs["spim"].path,
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
    resources:
        mem_mb=32000,
        runtime=lambda wildcards: int(200.0 / (3.0 ** float(wildcards.level))),  #rough estimate
    script:
        "../scripts/vesselfm.py"
