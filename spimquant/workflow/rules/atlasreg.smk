
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
        None
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
    shell:
        "c3d {input.corrected} {input.mask} -multiply -o {output.masked}"


rule affine_reg:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=config["atlasreg"]["stain"],
            level=config["atlasreg"]["level"],
            desc=config["atlasreg"]["desc"],
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
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -a -dof 12 -ia-image-centers -m NMI -o {output.xfm_ras} && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.xfm_ras}"


rule deform_reg:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=config["atlasreg"]["stain"],
            level=config["atlasreg"]["level"],
            desc=config["atlasreg"]["desc"],
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
        dseg=rules.import_dseg.output.dseg,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        label_tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
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
        chunks=(50, 50, 50),
        zooms=None,
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
        chunks=(50, 50, 50),
        zooms=None,
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
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        chunks=(50, 50, 50),
        zooms=None,
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
    script:
        "../scripts/deform_to_template_nii.py"


rule deform_to_template_nii_zoomed:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        chunks=(50, 50, 50),
        zooms=lambda wildcards: (
            float(wildcards.res) / 1000,
            float(wildcards.res) / 1000,
            float(wildcards.res) / 1000,
        ),  #None #same resolution as template if NOne
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
    threads: 32
    script:
        "../scripts/deform_to_template_nii.py"


rule deform_to_template_nii_nb:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        chunks=(20, 20, 20),
        zooms=None,  #same resolution as template if NOne
    output:
        nii=bids(
            root=root,
            datatype="micr",
            desc="deformnb",
            space="{template}",
            stain="{stain}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    container:
        None
    threads: 32
    notebook:
        "../notebooks/deform_to_template_nii.py.ipynb"


rule deform_transform_labels_to_subj:
    input:
        ref_ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        invwarp_nii=rules.deform_reg.output.invwarp,
        flo_nii=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
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
    script:
        "../scripts/deform_transform_channel_to_template_nii.py"


rule transform_labels_to_zoomed_template:
    input:
        dseg=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
        ref=bids(
            root=root,
            datatype="micr",
            desc="deform",
            space="{template}",
            stain=config["atlasreg"]["stain"],
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
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n NearestNeighbor "
        " -i {input.dseg} -o {output.dseg} "
        " -r {input.ref} "
