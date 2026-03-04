rule run_vesselfm:
    """Run Vesselfm
    """
    input:
        spim=inputs["spim"].path
    output:
        new_res=bids(
            root=root,
            datatype="new_res",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.ome.zarr",
            **inputs["spim"].wildcards,
        )
    conda:
        "../envs/vesselfm.yaml"
    shell:
        "python -m vesselfm.cli --input-folder {input.spim} --output-folder {output.new_res} --enable-dask-chunking --chunk-size 128 128 128 --downsample-level 0"