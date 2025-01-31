
rule n4:
    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),
    output:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4",
            suffix="biasfield.nii",
            **inputs["spim"].wildcards
        ),
    container:
        config["containers"]["ants"]
    shell:
        "N4BiasFieldCorrection -i {input.nii}"
        " -o [{output.corrected},{output.biasfield}]"
        " -x {input.mask} "
        " -d 3 -v "


rule apply_mask_to_corrected:
    input:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),
    output:
        masked=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="N4brain",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    container:
        config["containers"]["itksnap"]
    shell:
        "c3d {input.corrected} {input.mask} -multiply -o {output.masked}"


rule affine_reg:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    output:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="{template}",
            stain=stain_for_reg,
            desc="affinewarped",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    log:
        bids(
            root="logs",
            datatype="affine_reg",
            space="{template}",
            suffix="log.txt",
            **inputs["spim"].wildcards
        ),
    threads: 32
    container:
        config["containers"]["itksnap"]
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -a -dof 12 -ia-image-centers -m NMI -o {output.xfm_ras} && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.xfm_ras}"


rule convert_ras_to_itk:
    input:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards
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
            **inputs["spim"].wildcards
        ),
    container:
        config["containers"]["itksnap"]
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


rule deform_reg:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        xfm_ras=rules.affine_reg.output.xfm_ras,
    params:
        iters="100x50",
        metric="NMI",
        sigma1="4vox",
        sigma2="2vox",
    output:
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii",
            **inputs["spim"].wildcards
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            from_="{template}",
            to="subject",
            suffix="warp.nii",
            **inputs["spim"].wildcards
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="{template}",
            stain=stain_for_reg,
            desc="deformwarped",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    log:
        bids(
            root="logs",
            datatype="deform_reg",
            space="{template}",
            suffix="log.txt",
            **inputs["spim"].wildcards
        ),
    threads: 32
    resources:
        mem_mb=16000,
    container:
        config["containers"]["itksnap"]
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
                    **inputs["spim"].wildcards
                )
            )
        ),
    threads: 10
    log:
        bids(
            root="logs",
            datatype="resample_labels_to_zarr",
            space="{template}",
            suffix="log.txt",
            **inputs["spim"].wildcards
        ),
    script:
        "../scripts/resample_labels_to_zarr.py"


rule affine_zarr_to_template_nii:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        ref_opts={"chunks": (1, 50, 50, 50)},
    output:
        nii=bids(
            root=root,
            datatype="micr",
            desc="affine",
            space="{template}",
            stain="{stain}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    threads: 32
    script:
        "../scripts/affine_to_template_nii.py"


rule affine_zarr_to_template_ome_zarr:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        ref_opts={"chunks": (1, 50, 50, 50)},
    output:
        ome_zarr=directory(
            bids(
                root=root,
                datatype="micr",
                desc="affine",
                space="{template}",
                stain="{stain}",
                suffix="spim.ome.zarr",
                **inputs["spim"].wildcards
            )
        ),
    threads: 32
    script:
        "../scripts/affine_to_template_ome_zarr.py"


rule deform_zarr_to_template_nii:
    input:
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
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
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    threads: 32
    container:
        None
    script:
        "../scripts/deform_to_template_nii.py"


rule deform_to_template_nii_zoomed:
    input:
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        **get_storage_creds(inputs["spim"].path),
    params:
        ome_zarr=inputs["spim"].path,
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
    container: None
    output:
        nii=bids(
            root=root,
            datatype="micr",
            desc="deform",
            space="{template}",
            stain="{stain}",
            res="{res}um",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    resources: 
        mem_mb=15000
    threads: 4
    script:
        "../scripts/deform_to_template_nii.py"


rule deform_spim_nii_to_template_nii:
    input:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        ref=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"),
        xfm_itk=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards
        ),
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii",
            **inputs["spim"].wildcards
        ),
    output:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            space="{template}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    threads: 32
    container:
        config["containers"]["ants"]
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.spim} -o {output.spim} "
        " -r {input.ref} -t {input.warp} {input.xfm_itk}"


rule deform_template_dseg_to_subject_nii:
    """ use this to interpolate labels for each blob, and to calculate volumes"""
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
        xfm_ras=rules.init_affine_reg.output.xfm_ras,
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
            **inputs["spim"].wildcards
        ),
    threads: 32
    container:
        config["containers"]["itksnap"]
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        " -ri NN "
        "  -rm {input.dseg} {output.dseg} "
        "  -r {input.xfm_ras},-1 {input.invwarp}"
        #note: LABEL interpolation not possible with >1000 labels


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
                root=root,
                datatype="micr",
                desc="deform",
                space="subject",
                suffix="dseg.zarr",
                **inputs["spim"].wildcards
            )
        ),
    container:
        None
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
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            space="{template}",
            res="{res}um",
            suffix="dseg.nii",
            **inputs["spim"].wildcards
        ),
    threads: 32
    container:
        config["containers"]["ants"]
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n NearestNeighbor "
        " -i {input.dseg} -o {output.dseg} "
        " -r {input.ref} "
