wildcard_constraints:
    level="[0-9]+",
    template="[a-zA-Z0-9]+",


rule get_downsampled_nii:
    input:
        zarr=inputs["spim"].path,
    output:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    threads: 32
    script:
        "../scripts/ome_zarr_to_nii.py"


rule import_anat:
    input:
        anat=lambda wildcards: format(config["atlases"][wildcards.template]["anat"]),
    output:
        anat=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    log:
        bids_tpl(
            root="logs",
            datatype="import_anat",
            template="{template}",
            suffix="log.txt",
        ),
    shell:
        "cp {input} {output}"


rule import_dseg:
    input:
        dseg=lambda wildcards: format(config["atlases"][wildcards.template]["dseg"]),
    output:
        dseg=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
    log:
        bids_tpl(
            root="logs",
            datatype="import_dseg",
            template="{template}",
            suffix="log.txt",
        ),
    shell:
        "cp {input} {output}"


rule import_labelmapper_lut:
    input:
        json=lambda wildcards: format(config["atlases"][wildcards.template]["lut"]),
    output:
        tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
    log:
        bids_tpl(
            root="logs", datatype="import_lut", template="{template}", suffix="log.txt"
        ),
    script:
        "../scripts/import_labelmapper_lut.py"


rule make_itksnap_lut:
    input:
        tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
    output:
        lut=bids_tpl(
            root=root, template="{template}", desc="itksnap", suffix="labels.txt"
        ),
    script:
        "../scripts/lut_bids_to_itksnap.py"


rule ome_zarr_to_zipstore:
    """ generic rule to process any ome.zarr from work """
    input:
        zarr=f"{work}/{{prefix}}.ome.zarr",
    output:
        zarr_zip=f"{root}/{{prefix}}.ome.zarr.zip",
    log:
        "logs/ome_zarr_to_zipstore/{prefix}.log",
    group:
        "preproc"
    shell:
        "7z a -mx0 -tzip {output.zarr_zip} {input.zarr}/. &> {log}"


rule lateralize_atlas_labels:
    """splits atlas label nii into left and right, adding an offset to right hemi"""
    input:
        dseg=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
    params:
        offset_rh=10000,
    output:
        dseg=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"),
    shell:
        "c3d {input.dseg} -as SEG -cmv -pop -pop -threshold 50% inf 1 0 -as MASK_RH "
        " -push SEG -times -as SEG_RH "
        " -push MASK_RH -replace 1 0 0 1 -as MASK_LH "
        " -push SEG -times -as SEG_LH "
        " -push SEG_RH -binarize -scale {params.offset_rh} -push SEG_RH -add "
        " -push SEG_LH -add -type int -o {output.dseg}"


rule lateralize_atlas_tsv:
    """splits atlas label tsv into left and right, adding an offset to right hemi"""
    input:
        tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
    output:
        tsv=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.tsv"),
    script:
        "../scripts/lateralize_atlas_tsv.py"
