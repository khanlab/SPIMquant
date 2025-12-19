
rule gaussian_biasfield:
    """simple bias field correction with gaussian"""
    input:
        spim=inputs["spim"].path,
    params:
        proc_level=5,
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        corrected=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="correctedgaussian",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
        biasfield=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="gaussian",
                    suffix="biasfield.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=15,
    script:
        "../scripts/gaussian_biasfield.py"


rule n4_biasfield:
    """N4 bias field correction with antspyx"""
    input:
        spim=inputs["spim"].path,
    params:
        proc_level=5,
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        corrected=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="correctedn4",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
        biasfield=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="n4",
                    suffix="biasfield.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=500000,
        disk_mb=2097152,
        runtime=60,
    script:
        "../scripts/n4_biasfield.py"


rule destripe:
    input:
        spim=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="SPIM.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        destripe_kwargs={},
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        destriped=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="{desc}destripe",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=500000,
        runtime=60,
    script:
        "../scripts/destripe.py"
