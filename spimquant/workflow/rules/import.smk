

wildcard_constraints:
    level="[0-9]+",
    template="[a-zA-Z0-9]+",


rule get_downsampled_nii:
    input:
        spim=inputs["spim"].path,
    params:
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    script:
        "../scripts/ome_zarr_to_nii.py"


rule import_template_anat:
    input:
        anat=lambda wildcards: storage(
            ancient(resources_path(config["templates"][wildcards.template]["anat"]))
        ),
    output:
        anat=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    log:
        bids_tpl(
            root="logs",
            datatype="import_anat",
            template="{template}",
            suffix="log.txt",
        ),
    script:
        "../scripts/copy_nii.py"


rule import_mask:
    input:
        mask=lambda wildcards: storage(
            ancient(resources_path(config["templates"][wildcards.template]["mask"]))
        ),
    output:
        mask=bids_tpl(
            root=root, template="{template}", desc="brain", suffix="mask.nii.gz"
        ),
    log:
        bids_tpl(
            root="logs",
            datatype="import_mask",
            template="{template}",
            suffix="log.txt",
        ),
    script:
        "../scripts/copy_nii.py"


rule generic_lut_bids_to_itksnap:
    input:
        tsv="{prefix}_dseg.tsv",
    output:
        lut="{prefix}_dseg.itksnap.txt",
    script:
        "../scripts/lut_bids_to_itksnap.py"


rule import_dseg:
    input:
        dseg=lambda wildcards: storage(
            ancient(
                resources_path(
                    config["templates"][wildcards.template]["atlases"][wildcards.seg][
                        "dseg"
                    ]
                )
            )
        ),
    output:
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
    script:
        "../scripts/copy_nii.py"


rule import_lut_tsv:
    input:
        tsv=lambda wildcards: ancient(
            resources_path(
                config["templates"][wildcards.template]["atlases"][wildcards.seg][
                    "tsv"
                ]
            )
        ),
    output:
        tsv=bids_tpl(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    shell:
        "cp {input} {output}"


rule import_DSURQE_tsv:
    input:
        csv=storage(config["templates"]["DSURQE"]["atlases"]["all"]["custom_csv"]),
    output:
        tsv=bids_tpl(root=root, template="DSURQE", seg="all", suffix="dseg.tsv"),
    script:
        "../scripts/import_DSURQE_dseg_tsv.py"
