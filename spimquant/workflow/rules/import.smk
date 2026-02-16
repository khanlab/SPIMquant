"""
Data import and format conversion workflow for SPIMquant.

This module handles importing reference data (templates, atlases, masks) and converting
SPIM data from OME-Zarr format to NIfTI format at various resolution levels.

Key components:
1. Template anatomy import (from resources directory or remote URLs)
2. Brain mask import
3. Atlas segmentation (dseg) and label table (TSV) import
4. OME-Zarr to NIfTI conversion at specified downsampling levels
5. Label lookup table conversion for visualization tools

All imported files are organized following BIDS conventions with template-specific
subdirectories (tpl-{template}/).
"""


wildcard_constraints:
    level="[0-9]+",
    template="[a-zA-Z0-9]+",


rule get_downsampled_nii:
    """Convert OME-Zarr to NIfTI at specified resolution level.
    
    Extracts a single resolution level from the multi-scale OME-Zarr data and
    converts it to NIfTI format. The level parameter determines which pyramid
    level to extract (0=highest resolution). Handles orientation conversion
    based on configuration.
    """
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
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/ome_zarr_to_nii.py"


localrules:
    import_template_anat,
    import_mask,
    generic_lut_bids_to_itksnap,
    import_dseg,
    import_lut_tsv,
    import_DSURQE_tsv,


rule import_template_anat:
    """Import template anatomical image.
    
    Copies or downloads the template reference image specified in the configuration.
    Supports both local paths (relative to resources/) and remote URLs.
    Templates can include ABAv3, gubra, MBMv3, turone, and MouseIn.
    """
    input:
        anat=lambda wildcards: storage(
            ancient(resources_path(config["templates"][wildcards.template]["anat"]))
        ),
    output:
        anat=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
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
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
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
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/lut_bids_to_itksnap.py"


rule import_dseg:
    """Import atlas segmentation (dseg) file.
    
    Copies or downloads the atlas parcellation file for the specified template
    and segmentation scheme. The dseg file contains discrete labels corresponding
    to anatomical regions defined in the companion TSV file.
    """
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
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
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
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    shell:
        "cp {input} {output}"


rule import_DSURQE_tsv:
    input:
        csv=storage(config["templates"]["DSURQE"]["atlases"]["all"]["custom_csv"]),
    output:
        tsv=bids_tpl(root=root, template="DSURQE", seg="all", suffix="dseg.tsv"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/import_DSURQE_dseg_tsv.py"
