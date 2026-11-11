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
            stain=stain_for_reg,
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
        "subj_reg"
    conda:
        "../envs/ants.yaml"
    threads: 1
    resources:
        mem_mb=scaled_mem_mb(2000, 22000),
        runtime=scaled_runtime(10, 30),
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
            stain=stain_for_reg,
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
        "subj_reg"
    conda:
        "../envs/c3d.yaml"
    threads: 1
    resources:
        mem_mb=scaled_mem_mb(2000, 14000),
        runtime=10,
    shell:
        "c3d {input.corrected} {input.mask} -multiply -o {output.masked}"


rule crop_template:
    """Crop template to hemisphere for registration.

    Allows registration to left or right hemisphere only by cropping the template
    to the specified hemisphere region.
    """
    input:
        template=bids(root=root, template="{template}", suffix="{suffix}.nii.gz"),
    output:
        cropped=bids(
            root=root,
            template="{template}",
            desc="{hemisphere}crop",
            suffix="{suffix}.nii.gz",
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    params:
        hemisphere="{hemisphere}",
    script:
        "../scripts/crop_template.py"


def _reg_moving_stain(stain):
    return bids(
        root=root,
        datatype="micr",
        stain=stain,
        level=config["registration_level"],
        desc=config["templatereg"]["desc"],
        suffix="SPIM.nii.gz",
        **inputs["spim"].wildcards,
    )


_reg_moving = _reg_moving_stain(stain_for_reg)


def _reg_other_stains(wildcards):
    """Stains registered as well as stain_for_reg, for multichannel methods."""
    return reg_stains[1:] if reg_presets[wildcards.reg].get("multichannel") else []


def _reg_extra_fixed(wildcards):
    return [
        get_template_path(root, wildcards.template, config.get("template_crop"), stain)
        for stain in _reg_other_stains(wildcards)
    ]


def _reg_extra_moving(wildcards):
    return [
        _reg_moving_stain(stain).format(**wildcards)
        for stain in _reg_other_stains(wildcards)
    ]


def _reg_extra_pairs(wildcards):
    """greedy image pairs for the other stains, template and subject."""
    return " ".join(
        f"-i {fixed} {moving}"
        for fixed, moving in zip(
            _reg_extra_fixed(wildcards), _reg_extra_moving(wildcards)
        )
    )


def _reg_xfm(suffix, from_="subject", to="{template}", desc="{reg}", **kwargs):
    return bids(
        root=root,
        datatype="xfm",
        from_=from_,
        to=to,
        desc=desc,
        suffix=suffix,
        **kwargs,
        **inputs["spim"].wildcards,
    )


def _reg_warped(desc):
    return bids(
        root=root,
        datatype="xfm",
        space="{template}",
        stain=stain_for_reg,
        desc=desc,
        suffix="SPIM.nii.gz",
        **inputs["spim"].wildcards,
    )


rule affine_reg:
    """Affine registration with greedy, for methods with backend: greedy.

    Writes the transformation and the warped image for QC. Methods with
    affine: false use the identity instead. The same rule, with other paths,
    registers the original template to a cohort template (cohort_template.smk).
    """
    input:
        fixed=get_template_for_reg,
        moving=_reg_moving,
        extra_fixed=_reg_extra_fixed,
        extra_moving=_reg_extra_moving,
    output:
        xfm_ras=temp(_reg_xfm("xfm.txt", desc="{reg}affine", type_="ras")),
        warped=temp(_reg_warped("{reg}affinewarped")),
    log:
        bids(
            root="logs",
            datatype="affine_reg",
            space="{template}",
            desc="{reg}",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    wildcard_constraints:
        reg=reg_constraint("greedy"),
    group:
        "subj_reg"
    threads: 32
    resources:
        mem_mb=reg_mem_mb(16000),
        runtime=15,
    params:
        affine=lambda wildcards: reg_presets[wildcards.reg]["affine"],
        init="-ia-image-centers",
        metric=lambda wildcards: reg_presets[wildcards.reg]["metric"],
        iters=lambda wildcards: reg_iters(reg_presets[wildcards.reg], "affine_iters"),
        identity=reg_identity_ras,
        extra_pairs=_reg_extra_pairs,
    shell:
        "(if [ '{params.affine}' = 'True' ]; then "
        " greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} {params.extra_pairs} "
        "  -a -dof 12 {params.init} -m {params.metric} -o {output.xfm_ras} -n {params.iters}; "
        "else {params.identity} > {output.xfm_ras}; fi && "
        " greedy -threads {threads} -d 3 -rf {input.fixed} "
        "  -rm {input.moving} {output.warped} "
        "  -r {output.xfm_ras}) &> {log}"


rule convert_ras_to_itk:
    """Convert greedy's RAS affine to ITK format, for applying with ANTs tools."""
    input:
        xfm_ras=rules.affine_reg.output.xfm_ras,
    output:
        xfm_itk=temp(_reg_xfm("xfm.txt", desc="{reg}affine", type_="itk")),
    wildcard_constraints:
        reg=reg_constraint("greedy"),
    group:
        "subj_reg"
    conda:
        "../envs/c3d.yaml"
    threads: 1
    resources:
        mem_mb=1500,
        runtime=10,
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


rule deform_reg:
    """Deformable registration with greedy, starting from affine_reg.

    Writes the forward and inverse warps, and the warped image for QC.
    """
    input:
        fixed=get_template_for_reg,
        moving=_reg_moving,
        extra_fixed=_reg_extra_fixed,
        extra_moving=_reg_extra_moving,
        xfm_ras=rules.affine_reg.output.xfm_ras,
    output:
        warp=temp(_reg_xfm("warp.nii.gz")),
        invwarp=temp(_reg_xfm("warp.nii.gz", from_="{template}", to="subject")),
        warped=temp(_reg_warped("{reg}deformwarped")),
    log:
        bids(
            root="logs",
            datatype="deform_reg",
            space="{template}",
            desc="{reg}",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    wildcard_constraints:
        reg=reg_constraint("greedy"),
    group:
        "subj_reg"
    threads: 32
    resources:
        mem_mb=reg_mem_mb(20000),
        runtime=5 if config["sloppy"] else 30,
    params:
        iters=lambda wildcards: reg_iters(reg_presets[wildcards.reg], "deform_iters"),
        metric=lambda wildcards: reg_presets[wildcards.reg]["metric"],
        sigma1=lambda wildcards: reg_presets[wildcards.reg]["sigma1"],
        sigma2=lambda wildcards: reg_presets[wildcards.reg]["sigma2"],
        extra_pairs=_reg_extra_pairs,
    shell:
        "(greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} {params.extra_pairs} "
        " -it {input.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} && "
        " greedy -threads {threads} -d 3 -rf {input.fixed} "
        "  -rm {input.moving} {output.warped} "
        "  -r {output.warp} {input.xfm_ras}) &> {log}"


rule reg_ants:
    """Affine and deformable registration with antsRegistration (backend: ants)."""
    input:
        fixed=get_template_for_reg,
        moving=_reg_moving,
        extra_fixed=_reg_extra_fixed,
        extra_moving=_reg_extra_moving,
    output:
        xfm_itk=temp(_reg_xfm("xfm.txt", desc="{reg}affine", type_="itk")),
        warp=temp(_reg_xfm("warp.nii.gz")),
        invwarp=temp(_reg_xfm("warp.nii.gz", from_="{template}", to="subject")),
        warped_affine=temp(_reg_warped("{reg}affinewarped")),
        warped=temp(_reg_warped("{reg}deformwarped")),
    log:
        bids(
            root="logs",
            datatype="reg_ants",
            space="{template}",
            desc="{reg}",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    wildcard_constraints:
        reg=reg_constraint("ants"),
    conda:
        "../envs/ants.yaml"
    threads: 32
    resources:
        mem_mb=reg_mem_mb(32000),
        runtime=per_attempt(5 if config["sloppy"] else 240),
    params:
        args=lambda wildcards, input: ants_reg_args(
            reg_presets[wildcards.reg],
            input.fixed,
            input.moving,
            center_init=True,
            extra_fixed=_reg_extra_fixed(wildcards),
            extra_moving=_reg_extra_moving(wildcards),
        ),
        identity=reg_identity_itk,
    shell:
        "(export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads}; "
        "tmp=$(mktemp -d) && "
        "antsRegistration {params.args} -o $tmp/reg_ && "
        "if [ -f $tmp/reg_0GenericAffine.mat ]; then "
        " ConvertTransformFile 3 $tmp/reg_0GenericAffine.mat {output.xfm_itk}; "
        "else {params.identity} > {output.xfm_itk}; fi && "
        "mv $tmp/reg_*[0-9]Warp.nii.gz {output.warp} && "
        "mv $tmp/reg_*InverseWarp.nii.gz {output.invwarp} && "
        "antsApplyTransforms -d 3 -u float -i {input.moving} -r {input.fixed} "
        " -o {output.warped_affine} -t {output.xfm_itk} && "
        "antsApplyTransforms -d 3 -u float -i {input.moving} -r {input.fixed} "
        " -o {output.warped} -t {output.warp} {output.xfm_itk} && "
        "rm -r $tmp) &> {log}"


rule reg_fireants:
    """Affine and deformable registration with FireANTs (backend: fireants, GPU)."""
    input:
        fixed=get_template_for_reg,
        moving=_reg_moving,
        extra_fixed=_reg_extra_fixed,
        extra_moving=_reg_extra_moving,
    output:
        xfm_itk=temp(_reg_xfm("xfm.txt", desc="{reg}affine", type_="itk")),
        warp=temp(_reg_xfm("warp.nii.gz")),
        invwarp=temp(_reg_xfm("warp.nii.gz", from_="{template}", to="subject")),
        warped_affine=temp(_reg_warped("{reg}affinewarped")),
        warped=temp(_reg_warped("{reg}deformwarped")),
    wildcard_constraints:
        reg=reg_constraint("fireants"),
    threads: 8
    resources:
        gpu=1,
        tasks_per_gpu=0,
        cpus_per_gpu=8,
        mem_mb=reg_mem_mb(32000),
        runtime=30,
    params:
        preset=lambda wildcards: {
            **reg_presets[wildcards.reg],
            "affine_iters": reg_iters(reg_presets[wildcards.reg], "affine_iters"),
            "deform_iters": reg_iters(reg_presets[wildcards.reg], "deform_iters"),
        },
        init_centers=lambda wildcards: reg_presets[wildcards.reg]["affine"],
        extra_fixed=_reg_extra_fixed,
        extra_moving=_reg_extra_moving,
    script:
        "../scripts/reg_fireants.py"


rule compose_subject_to_template_warp:
    """Compose a method's affine and warp into a composite warp field per direction.

    from-subject_to-{template} resamples subject images into template space,
    and from-{template}_to-subject is its inverse.
    """
    input:
        ref=rules.import_template_anat.output.anat,
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        xfm_itk=_reg_xfm("xfm.txt", desc="{reg}affine", type_="itk"),
        warp=_reg_xfm("warp.nii.gz"),
        invwarp=_reg_xfm("warp.nii.gz", from_="{template}", to="subject"),
    output:
        xfm_composite=_reg_xfm("xfm.nii.gz"),
        xfm_composite_inv=_reg_xfm("xfm.nii.gz", from_="{template}", to="subject"),
    wildcard_constraints:
        reg=reg_constraint(),
    group:
        "subj_reg"
    conda:
        "../envs/ants.yaml"
    threads: 4
    resources:
        mem_mb=10000,
        runtime=15,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite},1] "
        " -r {input.ref} -t {input.warp} {input.xfm_itk} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite_inv},1] "
        " -r {input.subject} -t [{input.xfm_itk},1] {input.invwarp}"


rule compose_subject_warp_from_cohort_build:
    """Subject<->cohort template warps from the warps saved by the cohort build.

    Chains the subject's affine to the original template with its warp to the
    cohort template, instead of registering the subject again. Used when
    --template is a cohort template built from this subject. A point in the output
    space goes through the first transform listed after -t first.
    """
    input:
        ref=rules.import_template_anat.output.anat,
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        warp=lambda wildcards: _cohort_warp_paths(wildcards)["warp"],
        invwarp=lambda wildcards: _cohort_warp_paths(wildcards)["invwarp"],
        to_orig=lambda wildcards: _cohort_warp_paths(wildcards)["to_orig"],
        from_orig=lambda wildcards: _cohort_warp_paths(wildcards)["from_orig"],
    output:
        xfm_composite=_reg_xfm("xfm.nii.gz", desc="cohort"),
        xfm_composite_inv=_reg_xfm(
            "xfm.nii.gz", from_="{template}", to="subject", desc="cohort"
        ),
    group:
        "subj_reg"
    conda:
        "../envs/ants.yaml"
    threads: 4
    resources:
        mem_mb=10000,
        runtime=15,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite},1] "
        " -r {input.ref} -t {input.warp} {input.to_orig} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite_inv},1] "
        " -r {input.subject} -t {input.from_orig} {input.invwarp}"


rule select_subject_to_template_warp:
    """Link the warps every downstream rule uses to the subject's selected method."""
    input:
        xfm_composite=lambda wildcards: _reg_xfm("xfm.nii.gz").format(
            reg=subject_reg_method(wildcards), **wildcards
        ),
        xfm_composite_inv=lambda wildcards: _reg_xfm(
            "xfm.nii.gz", from_="{template}", to="subject"
        ).format(reg=subject_reg_method(wildcards), **wildcards),
    output:
        xfm_composite=bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to="{template}",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
        xfm_composite_inv=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to="subject",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
    localrule: True
    shell:
        "ln -srf {input.xfm_composite} {output.xfm_composite} && "
        "ln -srf {input.xfm_composite_inv} {output.xfm_composite_inv}"


rule deform_zarr_to_template_nii:
    input:
        xfm_composite=bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to="{template}",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
        ref_nii=get_template_for_reg,
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
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        ome_zarr=spim_input,
        flo_opts={"level": 2},  #downsampling level to use (TODO: set this automatically based on ref resolution?)
        do_downsample=True,  #whether to perform further downsampling before transforming
        downsample_opts={"along_z": 4},  #could also be determined automatically 
        ref_opts={"chunks": (1, 100, 100, 100)},
    script:
        "../scripts/deform_to_template_nii.py"


rule deform_to_template_nii_zoomed:
    input:
        ome_zarr=spim_input,
        spim_json=spim_json_input,
        xfm_composite=bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to="{template}",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
        ref_nii=get_template_for_reg,
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
    threads: 4
    resources:
        mem_mb=16000,
        runtime=15,
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
        xfm_composite=bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to="{template}",
            suffix="xfm.nii.gz",
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
        "subj_reg"
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=scaled_mem_mb(8000, 32000),
        runtime=10,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.spim} -o {output.spim} "
        " -r {input.ref} -t {input.xfm_composite}"


rule deform_template_dseg_to_subject_nii:
    """Transform template atlas labels to subject space.

    Applies composite inverse warp to bring template segmentation labels into
    subject space. Uses nearest-neighbor interpolation to preserve discrete label
    values. This enables atlas-based analysis in native subject space.
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
        dseg=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"),
        xfm_composite_inv=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to="subject",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        dseg=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level="{level}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_reg"
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=scaled_mem_mb(2500, 16000),
        runtime=10,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n NearestNeighbor "
        " -i {input.dseg} -o {output.dseg} "
        " -r {input.ref} -t {input.xfm_composite_inv}"


_qc_warped = bids(
    root=root,
    datatype="xfm",
    space="{template}",
    stain="{stain}",
    desc="{reg}qcwarped",
    suffix="SPIM.nii.gz",
    **inputs["spim"].wildcards,
)
# the template stains in the reports, or else the one compared to the template's anat
_report_stains = reg_stains or [stain_for_reg]


def _report_fixed(wildcards, stain):
    suffix = stain if use_spim_template else None
    return get_template_path(
        root, wildcards.template, config.get("template_crop"), suffix
    )


rule warp_stain_for_report:
    """Resample a stain with a method's composite warp, for the reports."""
    input:
        moving=lambda wildcards: _reg_moving_stain(wildcards.stain).format(**wildcards),
        ref=lambda wildcards: _report_fixed(wildcards, wildcards.stain),
        xfm=_reg_xfm("xfm.nii.gz"),
    output:
        temp(_qc_warped),
    wildcard_constraints:
        reg=reg_constraint(),
        stain="|".join(_report_stains),
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=per_attempt(16000),
        runtime=per_attempt(15),
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.moving} -o {output} -r {input.ref} -t {input.xfm}"


rule registration_qc_report:
    """Report on the affine and deformable stages of one registration method."""
    input:
        template=get_template_for_reg,
        warped_affine=_reg_warped("{reg}affinewarped"),
        warped_deform=_reg_warped("{reg}deformwarped"),
        warp=_reg_xfm("xfm.nii.gz"),
        extra_fixed=lambda wildcards: [
            _report_fixed(wildcards, stain) for stain in _report_stains[1:]
        ],
        extra_warped=lambda wildcards: [
            _qc_warped.format(stain=stain, **wildcards) for stain in _report_stains[1:]
        ],
        dseg=lambda wildcards: bids(
            root=root,
            template=wildcards.template,
            seg=list(atlas_segs)[0],
            suffix="dseg.nii.gz",
        ),
    output:
        report_html=_reg_xfm("regqc.html", from_="SPIM", stain=stain_for_reg),
    wildcard_constraints:
        reg=reg_constraint(),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        subject=lambda wildcards: f"sub-{wildcards.subject}",
        template=lambda wildcards: f"tpl-{wildcards.template}",
        stains=_report_stains,
        method=lambda wildcards: wildcards.reg,
        parameters=lambda wildcards: reg_note(reg_presets[wildcards.reg]),
    script:
        "../scripts/reg_qc_report.py"


rule reg_tune_report:
    """Compare registration methods for one subject, using the same report as QC."""
    input:
        fixed=get_template_for_reg,
        extra_fixed=lambda wildcards: [
            _report_fixed(wildcards, stain) for stain in _report_stains[1:]
        ],
        # from the composite warps, which are kept, so that the report can be
        # written again without repeating the registrations
        warped=lambda wildcards: [
            _qc_warped.format(reg=reg, stain=stain, **wildcards)
            for stain in _report_stains
            for reg in config["registration_tune_methods"]
        ],
        warps=lambda wildcards: [
            _reg_xfm("xfm.nii.gz").format(reg=reg, **wildcards)
            for reg in config["registration_tune_methods"]
        ],
        labels=lambda wildcards: [
            bids(
                root=root,
                template=wildcards.template,
                seg=list(atlas_segs)[0],
                suffix="dseg.nii.gz",
            )
        ]
        * len(config["registration_tune_methods"]),
    output:
        report_html=_reg_xfm(
            "regqc.html", from_="SPIM", stain=stain_for_reg, desc="tune"
        ),
    threads: 1
    resources:
        mem_mb=32000,
        runtime=60,
    params:
        subject=lambda wildcards: f"sub-{wildcards.subject}",
        template=lambda wildcards: f"tpl-{wildcards.template}",
        stains=_report_stains,
        methods=config["registration_tune_methods"],
        selected=config["registration_method"],
        notes={m: reg_note(reg_presets[m]) for m in reg_presets},
    script:
        "../scripts/reg_tune_report.py"


rule all_reg_tune:
    """Compare registration methods per subject (--registration_tune_methods).

    Invoke with: ``spimquant <bids_dir> <out_dir> tune_reg``
    """
    input:
        ref_spim.expand(
            rules.reg_tune_report.output.report_html,
            template=config["template"],
        ),
