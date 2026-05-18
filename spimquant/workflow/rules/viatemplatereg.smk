"""
Via-template registration workflow for SPIMquant.

This module handles multi-step template registration where subjects are first
registered to an intermediate template (via_template), and then the intermediate
template is registered to the final template. The transforms are concatenated
to produce a single composite warp from subject to final template space.

This is enabled with the --via-template option and is useful when direct
registration from subject SPIM to a target template is challenging due to
large appearance differences, by using an intermediate template that is more
similar to the SPIM data.

Key workflow stages:
1. Affine registration from via-template anatomy to final template anatomy
2. Convert affine to ITK format
3. Deformable registration from via-template to final template
4. Compose via-template → final template transforms into a single composite
5. Compose subject → via-template and via-template → final template composites
   into a single subject → final template composite warp
"""


rule affine_reg_template_to_template:
    """Perform affine registration from via-template to final template using greedy.

    Registers the via-template anatomy to the final template anatomy using 12-DOF
    affine transformation. Used as the first step in template-to-template registration.
    """
    input:
        template=rules.import_template_anat.output.anat,
        via_template_anat=bids(root=root, template=via_template, suffix="anat.nii.gz"),
    params:
        iters="10x0x0" if config["sloppy"] else "100x100",
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="xfm",
                from_=via_template,
                to="{template}",
                type_="ras",
                desc="affine",
                suffix="xfm.txt",
            )
        ),
        warped=temp(
            bids(
                root=root,
                datatype="xfm",
                space="{template}",
                from_=via_template,
                desc="affinewarped",
                suffix="anat.nii.gz",
            )
        ),
    log:
        bids(
            root="logs",
            datatype="affine_reg_t2t",
            space="{template}",
            suffix="log.txt",
        ),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.via_template_anat} "
        " -a -dof 12 -ia-image-centers -m NMI -o {output.xfm_ras} -n {params.iters} "
        " > {log} 2>&1 && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.via_template_anat} {output.warped} "
        "  -r {output.xfm_ras} >> {log} 2>&1"


rule convert_ras_to_itk_template_to_template:
    """Convert RAS affine transform (via-template to template) to ITK format.

    Converts greedy's RAS affine format to ITK format for compatibility with
    ANTs tools used in the subsequent compose step.
    """
    input:
        xfm_ras=bids(
            root=root,
            datatype="xfm",
            from_=via_template,
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
        ),
    output:
        xfm_itk=temp(
            bids(
                root=root,
                datatype="xfm",
                from_=via_template,
                to="{template}",
                type_="itk",
                desc="affine",
                suffix="xfm.txt",
            )
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


rule deform_reg_template_to_template:
    """Perform deformable registration from via-template to final template using greedy.

    Refines the affine registration with deformable (non-linear) registration.
    Outputs forward warp, inverse warp, and warped image for QC.
    """
    input:
        template=rules.import_template_anat.output.anat,
        via_template_anat=bids(root=root, template=via_template, suffix="anat.nii.gz"),
        xfm_ras=bids(
            root=root,
            datatype="xfm",
            from_=via_template,
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
        ),
    params:
        iters="10x0x0" if config["sloppy"] else "100x50",
        metric="NMI",
        sigma1="4vox",
        sigma2="2vox",
    output:
        warp=temp(
            bids(
                root=root,
                datatype="xfm",
                from_=via_template,
                to="{template}",
                suffix="warp.nii.gz",
            )
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="xfm",
                from_="{template}",
                to=via_template,
                suffix="warp.nii.gz",
            )
        ),
        warped=temp(
            bids(
                root=root,
                datatype="xfm",
                space="{template}",
                from_=via_template,
                desc="deformwarped",
                suffix="anat.nii.gz",
            )
        ),
    log:
        bids(
            root="logs",
            datatype="deform_reg_t2t",
            space="{template}",
            suffix="log.txt",
        ),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5 if config["sloppy"] else 30,
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.via_template_anat} "
        " -it {input.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} "
        " > {log} 2>&1 && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.via_template_anat} {output.warped} "
        "  -r {output.warp} {input.xfm_ras} >> {log} 2>&1"


rule compose_template_to_template_warp:
    """Compose affine and deformable transforms for via-template to final template.

    Creates composite transformation fields by concatenating the affine and
    deformable registration transforms for the template-to-template registration.
    This simplifies downstream applications that combine the subject-to-via and
    via-to-template transforms.
    """
    input:
        ref=rules.import_template_anat.output.anat,
        via_template_anat=bids(root=root, template=via_template, suffix="anat.nii.gz"),
        xfm_itk=bids(
            root=root,
            datatype="xfm",
            from_=via_template,
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
        ),
        warp=bids(
            root=root,
            datatype="xfm",
            from_=via_template,
            to="{template}",
            suffix="warp.nii.gz",
        ),
        invwarp=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to=via_template,
            suffix="warp.nii.gz",
        ),
    output:
        xfm_composite=bids(
            root=root,
            datatype="xfm",
            from_=via_template,
            to="{template}",
            suffix="xfm.nii.gz",
        ),
        xfm_composite_inv=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to=via_template,
            suffix="xfm.nii.gz",
        ),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite},1] "
        " -r {input.ref} -t {input.warp} {input.xfm_itk} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite_inv},1] "
        " -r {input.via_template_anat} -t {input.invwarp} [{input.xfm_itk},1]"


rule compose_subject_via_to_template_warp:
    """Compose subject-to-via and via-to-template warps into subject-to-template.

    Creates the composite transformation from subject space to the final template
    space by concatenating:
    1. subject → via_template (from compose_subject_to_template_warp with template=via_template)
    2. via_template → template (from compose_template_to_template_warp)

    Also creates the inverse composite (template → subject) by concatenating the
    corresponding inverse transforms in reverse order.
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
        xfm_composite_subj_to_via=bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to=via_template,
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
        xfm_composite_via_to_template=bids(
            root=root,
            datatype="xfm",
            from_=via_template,
            to="{template}",
            suffix="xfm.nii.gz",
        ),
        xfm_composite_inv_via_to_subj=bids(
            root=root,
            datatype="xfm",
            from_=via_template,
            to="subject",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
        xfm_composite_inv_template_to_via=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to=via_template,
            suffix="xfm.nii.gz",
        ),
    output:
        xfm_composite=bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to="{template}",
            via=via_template,
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
        xfm_composite_inv=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to="subject",
            via=via_template,
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite},1] "
        " -r {input.ref} "
        " -t {input.xfm_composite_via_to_template} "
        " -t {input.xfm_composite_subj_to_via} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite_inv},1] "
        " -r {input.subject} "
        " -t {input.xfm_composite_inv_via_to_subj} "
        " -t {input.xfm_composite_inv_template_to_via}"
