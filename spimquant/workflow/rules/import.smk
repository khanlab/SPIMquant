wildcard_constraints:
    level="[0-9]+",
    template="[a-zA-Z0-9]+",


rule get_downsampled_nii:
    """Downsample the (ch, z, y, x) SPIM OME zarr file, and save the result to a nifti file.
     
    Downsampling is done over x, y axis by taking the image under 'image.ome.zarr/level', and the along z 
    further to get a downsampled .nii image.
    
    input:
        zarr: The input brain scan OME zarr image (float32)
    
    output:
        nii: The down-sampled brain scan image in .nii format (float64)
    """
    input:
        zarr=cconfig.inputs["spim"].path,
    output:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii",
            **cconfig.inputs["spim"].wildcards
        ),
    threads: 32
    script:
        "../scripts/ome_zarr_to_nii.py"


rule import_anat:
    """Copy the anat .nii file into the `root` directory
    
    Usually, this is the image available (uploaded to github) in spimquant/resources/ABAv3/P56_Atlas.nii.gz,
    (which you can view using the napari-nifti plugin) 
    
    input:
        dseg: A uint8 nifti file providing segmentation of the regions

    output:
        dseg: Same mask copied to root directory
    """
    input:
        anat=lambda wildcards: format(config["templates"][wildcards.template]["anat"]),
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
    """Copy the dseg .nii file into the `root` directory
    
    Usually, this is the image available (uploaded to github) in spimquant/resources/ABAv3/P56_Annotation.nii.gz 
    (which you can view using the napari-nifti plugin)
    
    input:
        dseg: An int32 nifti file providing segmentation of the regions

    output:
        dseg: Same mask copied to root directory
    """
    input:
        dseg=lambda wildcards: format(config["templates"][wildcards.template]["dseg"]),
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
        json=lambda wildcards: format(config["templates"][wildcards.template]["lut"]),
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
        tsv=bids_tpl(root=root, template="{template}", desc="{desc}", suffix="dseg.tsv"),
    output:
        lut=bids_tpl(
            root=root, template="{template}", desc="{desc}", suffix="dseg.itksnap.txt"
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
    """splits atlas label nii into left and right, adding an offset to right hemi
    
    The output is a (z, y, x) instance segmentation mask of the same size as input, 
    but whose right (x positive) side objects are assigned numbers offset_rh higher number 
    than those on the left. e.g. An object labeled 800 on the left would be labeled 10800 
    on the right if offset_rh is 10000.
    
    input:
        dseg: Input mask is an int32 nifti file providing segmentation of the regions.

    output:
        dseg: A mask is of int32 type. The labels are the same as input but those on the right 
            have been incremented by 10000.
    """
    input:
        dseg=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
    params:
        offset_rh=10000,
    output:
        dseg=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"),
    container:
        config["containers"]["itksnap"]
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
