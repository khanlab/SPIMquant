rule cellpose:
    """ haven't tried on PI data or full-res, but so far doesn't really work well.."""
    input:
        zarr=inputs["spim"].path,
    params:
        level=lambda wildcards: int(wildcards.level),  #downsample-level to perform segmentation on
        chunks=(200, 200, 200),
    output:
        zarr=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="cellpose",
                suffix="dseg.zarr",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 6
    script:
        "../scripts/cellpose.py"
