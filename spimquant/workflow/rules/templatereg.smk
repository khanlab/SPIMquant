"""
Template registration workflow for SPIMquant.

This module handles the complete registration pipeline from subject SPIM data to template space.
It includes intensity correction, affine and deformable registration, transformation of labels,
and quality control reporting.

Key workflow stages:
1. Intensity correction (N4 bias field correction)
2. Brain masking and cropping
3. Affine registration to template
4. Deformable registration refinement
5. Forward/inverse transformations
6. Label resampling and warping
7. Quality control report generation

The workflow supports multiple templates (ABAv3, gubra, MBMv3, turone, MouseIn) and handles
both NIfTI and OME-Zarr formats for different resolution levels.
"""


rule n4:
    """Apply N4 bias field correction to SPIM images.
    
    Uses ANTs N4BiasFieldCorrection to correct intensity inhomogeneities within
    the brain mask. Outputs both the corrected image and the estimated bias field.
    """
    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4",
            suffix="biasfield.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    conda:
        "../envs/ants.yaml"
    shell:
        "N4BiasFieldCorrection -i {input.nii}"
        " -o [{output.corrected},{output.biasfield}]"
        " -x {input.mask} "
        " -d 3 -v "


rule apply_mask_to_corrected:
    """Apply brain mask to N4-corrected image.
    
    Multiplies the N4-corrected image with the brain mask to extract brain tissue,
    removing background and non-brain regions.
    """
    input:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        masked=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4brain",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.corrected} {input.mask} -multiply -o {output.masked}"


rule crop_template:
    """Crop template to hemisphere for registration.
    
    Allows registration to left or right hemisphere only by cropping the template
    to the specified hemisphere region.
    """
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    output:
        cropped=bids_tpl(
            root=root,
            template="{template}",
            desc="{hemisphere}crop",
            suffix="anat.nii.gz",
        ),
    params:
        hemisphere="{hemisphere}",
    threads: 1
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/crop_template.py"


rule affine_reg:
    """Perform affine registration to template using greedy.
    
    Registers subject SPIM to template space using 12-DOF affine transformation.
    Uses normalized mutual information (NMI) as the similarity metric and outputs
    both the transformation matrix and warped image for QC.
    """
    input:
        template=get_template_for_reg,
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    params:
        iters="10x0x0" if config["sloppy"] else "100x100",
    output:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="{template}",
            stain=stain_for_reg,
            desc="affinewarped",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    log:
        bids(
            root="logs",
            datatype="affine_reg",
            space="{template}",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -a -dof 12 -ia-image-centers -m NMI -o {output.xfm_ras} -n {params.iters} && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.xfm_ras}"


rule convert_ras_to_itk:
    """Convert RAS affine transform to ITK format.
    
    Converts greedy's RAS affine format to ITK format for compatibility with
    ANTs tools (used for applying transforms).
    """
    input:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
    output:
        xfm_itk=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


rule deform_reg:
    """Perform deformable registration to template using greedy.
    
    Refines affine registration with deformable (non-linear) registration.
    Outputs forward warp, inverse warp, and warped image. The deformation field
    captures local anatomical variations not handled by affine transformation.
    """
    input:
        template=get_template_for_reg,
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        xfm_ras=rules.affine_reg.output.xfm_ras,
    params:
        iters="10x0x0" if config["sloppy"] else "100x50",
        metric="NMI",
        sigma1="4vox",
        sigma2="2vox",
    output:
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii.gz",
            **inputs["spim"].wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            from_="{template}",
            to="subject",
            suffix="warp.nii.gz",
            **inputs["spim"].wildcards,
        ),
        warped=temp(
            bids(
                root=root,
                datatype="warps",
                space="{template}",
                stain=stain_for_reg,
                desc="deformwarped",
                suffix="SPIM.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj"
    log:
        bids(
            root="logs",
            datatype="deform_reg",
            space="{template}",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -it {input.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.warp} {input.xfm_ras} "


rule resample_labels_to_zarr:
    """TODO: add required OME metadata"""
    input:
        dseg=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"),
        label_tsv=bids_tpl(
            root=root, template="{template}", desc="LR", suffix="dseg.tsv"
        ),
        xfm_ras=rules.affine_reg.output.xfm_ras,
        zarr_zip=inputs["spim"].path,
    params:
        level_to_resample_to=0,
        max_downsampling_layers=config["ome_zarr"]["max_downsampling_layers"],
        label_name="dseg",
        scaling_method="nearest",
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    desc="resampled",
                    from_="{template}",
                    suffix="dseg.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 10
    resources:
        mem_mb=16000,
        disk_mb=2097152,
        runtime=15,
    log:
        bids(
            root="logs",
            datatype="resample_labels_to_zarr",
            space="{template}",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    script:
        "../scripts/resample_labels_to_zarr.py"


rule affine_zarr_to_template_nii:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        ref_nii=get_template_for_reg,
    params:
        ref_opts={"chunks": (1, 50, 50, 50)},
    output:
        nii=bids(
            root=root,
            datatype="micr",
            desc="affine",
            space="{template}",
            stain="{stain}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/affine_to_template_nii.py"


rule affine_zarr_to_template_ome_zarr:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        ref_nii=get_template_for_reg,
    params:
        ref_opts={"chunks": (1, 50, 50, 50)},
    output:
        ome_zarr=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    desc="affine",
                    space="{template}",
                    stain="{stain}",
                    suffix="spim.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        disk_mb=2097152,
        runtime=15,
    script:
        "../scripts/affine_to_template_ome_zarr.py"


rule deform_zarr_to_template_nii:
    input:
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=get_template_for_reg,
    params:
        ome_zarr=inputs["spim"].path,
        flo_opts={"level": 2},  #downsampling level to use (TODO: set this automatically based on ref resolution?)
        do_downsample=True,  #whether to perform further downsampling before transforming
        downsample_opts={"along_z": 4},  #could also be determined automatically 
        ref_opts={"chunks": (1, 100, 100, 100)},
    output:
        nii=bids(
            root=root,
            datatype="micr",
            desc="deform",
            space="{template}",
            stain="{stain}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/deform_to_template_nii.py"


rule deform_to_template_nii_zoomed:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=get_template_for_reg,
    params:
        flo_opts={},  #any additional flo znimg options
        do_downsample=True,  #whether to perform further downsampling before transforming
        downsample_opts={"along_z": 4},  #could also be determined automatically 
        ref_opts=lambda wildcards: {
            "chunks": (1, 50, 50, 50),
            "zooms": (
                float(wildcards.res) / 1000,
                float(wildcards.res) / 1000,
                float(wildcards.res) / 1000,
            ),
        },
    output:
        nii=bids(
            root=root,
            datatype="micr",
            desc="deform",
            space="{template}",
            stain="{stain}",
            res="{res}um",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 4
    resources:
        mem_mb=15000,
        runtime=15,
    script:
        "../scripts/deform_to_template_nii.py"


rule deform_spim_nii_to_template_nii:
    input:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        ref=rules.import_template_anat.output.anat,
        xfm_itk=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            space="{template}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.spim} -o {output.spim} "
        " -r {input.ref} -t {input.warp} {input.xfm_itk}"


rule deform_template_dseg_to_subject_nii:
    """Transform template atlas labels to subject space.
    
    Applies inverse warp to bring template segmentation labels into subject space.
    Uses nearest-neighbor interpolation to preserve discrete label values.
    This enables atlas-based analysis in native subject space.
    """
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
        xfm_itk=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
        invwarp=rules.deform_reg.output.invwarp,
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{level}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n NearestNeighbor "
        " -i {input.dseg} -o {output.dseg} "
        " -r {input.ref} -t  [{input.xfm_itk},1] {input.invwarp}"


rule copy_template_dseg_tsv:
    input:
        dseg=bids_tpl(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{level}",
            stain="{stain}",
            from_="{template}",
            suffix="dseg.tsv",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    shell:
        "cp {input} {output}"


""" this rule needs updating - use atlas/seg wildcard and proper script
rule deform_transform_labels_to_subj:
    input:
        ref_ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        invwarp_nii=rules.deform_reg.output.invwarp,
        flo_nii=bids_tpl(
            root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"
        ),
    output:
        zarr=directory(
            bids(
                root=work,
                datatype="micr",
                desc="deform",
                space="subject",
                suffix="dseg.zarr",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj"
    threads: 32
    script:  #TODO this script doesn't exist??
        "../scripts/deform_transform_channel_to_template_nii.py"

rule transform_labels_to_zoomed_template:
    input:
        dseg=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"),
        ref=bids(
            root=root,
            datatype="micr",
            desc="deform",
            space="{template}",
            stain=stain_for_reg,
            res="{res}um",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            space="{template}",
            res="{res}um",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n NearestNeighbor "
        " -i {input.dseg} -o {output.dseg} "
        " -r {input.ref} "
"""


rule registration_qc_report:
    """Generate registration quality control notebook with visualizations"""
    input:
        template=get_template_for_reg,
        subject=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        warped_affine=bids(
            root=root,
            datatype="warps",
            space="{template}",
            stain="{stain}",
            desc="affinewarped",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        warped_deform=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=config["registration_level"],
            space="{template}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii.gz",
            **inputs["spim"].wildcards,
        ),
        dseg=lambda wildcards: bids_tpl(
            root=root,
            template=wildcards.template,
            seg=list(atlas_segs)[0],
            suffix="dseg.nii.gz",
        ),
    output:
        report_html=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            space="{template}",
            suffix="regqc.html",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=10,
    script:
        "../scripts/reg_qc_report.py"
