rule run_vesselfm:
    input:
        spim=inputs["spim"].path,
    output:
        mask=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="vesselfm",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    script:
        "../scripts/vesselfm.py"
