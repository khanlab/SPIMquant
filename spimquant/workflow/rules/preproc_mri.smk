"""
MRI preprocessing and cross-modality registration workflow for SPIMquant.

This module handles co-registration of in-vivo MRI (T2w) with ex-vivo SPIM data,
enabling multi-modal analysis and assessment of tissue changes due to perfusion
fixation and optical clearing.

Key workflow stages:
1. N4 bias field correction of MRI
2. MRI to template registration (for brain masking)
3. Template brain mask to MRI transformation
4. MRI brain extraction
5. MRI to SPIM rigid+deformable registration
6. Parameter tuning rules for optimization
7. Concatenated transformations (MRI -> SPIM -> Template)
8. Jacobian determinant calculation (tissue deformation quantification)

This workflow is optional and used when both MRI and SPIM data are available
for the same subject. The MRI provides complementary anatomical information and
enables assessment of tissue shrinkage/expansion during sample preparation.
"""


def select_single_t2w(wildcards):

    files = inputs["T2w"].filter(subject=wildcards.subject).expand()
    if len(files) == 1:
        return files[0]
    elif len(files) == 0:
        raise ValueError(f"No T2w images found for f{wildcards}")
    else:
        raise ValueError(
            f"Multiple T2w images found for f{wildcards}, use --filter-T2w to select a single image"
        )


rule n4_mri:
    """Apply N4 bias field correction to T2w MRI.
    
    Uses ANTs N4BiasFieldCorrection to correct intensity inhomogeneities in
    the MRI image, improving subsequent registration performance.
    """
    input:
        nii=select_single_t2w,
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    shell:
        "N4BiasFieldCorrection -i {input.nii}"
        " -o {output.nii}"
        " -d 3 -v "


rule rigid_greedy_reg_mri_to_template:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
    params:
        iters="{iters}",  #"100x50x50",
        metric="NCC {radius}",  #3x3x3",
        sigma1="{gradsigma}vox",  #3
        sigma2="{warpsigma}vox",  #2
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="warps",
                from_="mri",
                to="{template}",
                type_="ras",
                desc="rigid",
                suffix="xfm.txt",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["T2w"].wildcards,
            )
        ),
        warp=temp(
            bids(
                root=root,
                datatype="warps",
                from_="mri",
                to="{template}",
                suffix="warp.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["T2w"].wildcards,
            )
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="warps",
                from_="{template}",
                to="mri",
                suffix="warp.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["T2w"].wildcards,
            )
        ),
        warped=temp(
            bids(
                root=root,
                datatype="warps",
                space="{template}",
                desc="deformwarped",
                suffix="T2w.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["T2w"].wildcards,
            )
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    shell:
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -a -dof 6 -ia-image-centers -m {params.metric} -o {output.xfm_ras} && "
        "greedy -threads {threads} -d 3 -i {input.template} {input.subject} "
        " -it {output.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} && "
        " greedy -threads {threads} -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.warp} {output.xfm_ras} "


#    log:
#        bids(
#            root="logs",
#            datatype="mri_rigid_greedy",
#            space="{template}",
#            suffix="log.txt",
#            **inputs["T2w"].wildcards,
#        ),


rule all_tune_mri_mask:
    input:
        inputs["T2w"].expand(
            bids(
                root=root,
                datatype="anat",
                desc="brain",
                suffix="mask.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["T2w"].wildcards,
            ),
            iters="100x50x50",
            gradsigma=range(3, 6),
            warpsigma=range(3, 6),
            radius=[f"{i}x{i}x{i}" for i in range(2, 5)],
        ),
    group:
        "subj"


rule transform_template_mask_to_mri:
    input:
        mask=bids_tpl(
            root=root,
            template=config["template_mri"],
            desc="brain",
            suffix="mask.nii.gz",
        ),
        ref=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_="mri",
            to=config["template_mri"],
            type_="ras",
            desc="rigid",
            suffix="xfm.txt",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs["T2w"].wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            from_=config["template_mri"],
            to="mri",
            suffix="warp.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs["T2w"].wildcards,
        ),
    output:
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs["T2w"].wildcards,
        ),
    shadow:
        "minimal"
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        " c3d_affine_tool {input.xfm_ras} -inv -o inv_rigid.txt && "
        " greedy -threads {threads} -d 3 -ri NN -rf {input.ref} "
        "  -rm {input.mask} {output.mask} "
        "  -r inv_rigid.txt {input.invwarp}"


rule apply_mri_brain_mask:
    input:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            iters="100x100x50",
            radius="2x2x2",
            gradsigma="3",
            warpsigma="3",
            **inputs["T2w"].wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4brain",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.nii} {input.mask} -multiply -resample 300% -o {output.nii}"


rule rigid_greedy_reg_mri_to_spim:
    """Register MRI to SPIM space using rigid + deformable registration.
    
    Performs initial rigid (6 or 12 DOF) alignment followed by deformable
    registration to align in-vivo MRI with ex-vivo SPIM data. Multiple parameter
    combinations can be tested via wildcards for optimization. Outputs both
    linear-only and deformable-warped results.
    """
    input:
        mri=bids(
            root=root,
            datatype="anat",
            desc="N4brain",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
        spim=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    params:
        iters="{iters}",  #"100x50x50",
        metric_rigid="NMI",
        dof="{dof}",
        metric="NCC {radius}",  #3x3x3",
        sigma1="{gradsigma}vox",  #3
        sigma2="{warpsigma}vox",  #2
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="warps",
                from_="mri",
                to="spim",
                type_="ras",
                desc="rigid",
                suffix="xfm.txt",
                iters="{iters}",
                dof="{dof}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["spim"].wildcards,
            )
        ),
        warp=temp(
            bids(
                root=root,
                datatype="warps",
                from_="mri",
                to="spim",
                suffix="warp.nii.gz",
                iters="{iters}",
                dof="{dof}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["spim"].wildcards,
            )
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="warps",
                from_="spim",
                to="mri",
                suffix="warp.nii.gz",
                iters="{iters}",
                dof="{dof}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["spim"].wildcards,
            )
        ),
        linwarped=bids(
            root=root,
            datatype="warps",
            space="SPIM",
            desc="linearwarped",
            suffix="T2w.nii.gz",
            iters="{iters}",
            dof="{dof}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs["spim"].wildcards,
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="SPIM",
            desc="deformwarped",
            suffix="T2w.nii.gz",
            iters="{iters}",
            dof="{dof}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    shell:
        "greedy -threads {threads} -d 3 -i {input.spim} {input.mri} "
        " -a -dof {params.dof} -ia-image-centers -m {params.metric_rigid} -o {output.xfm_ras} && "
        "greedy -threads {threads} -d 3 -i {input.spim} {input.mri} "
        " -it {output.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} && "
        " greedy -threads {threads} -d 3 -rf {input.spim} "
        "  -rm {input.mri} {output.warped} "
        "  -r {output.warp} {output.xfm_ras} && "
        " greedy -threads {threads} -d 3 -rf {input.spim} "
        "  -rm {input.mri} {output.linwarped} "
        "  -r {output.xfm_ras} "


#   log:
#       bids(
#           root="logs",
#           datatype="mri_spim_rigid_greedy",
#           suffix="log.txt",
#           **inputs["spim"].wildcards
#       ),


rule all_tune_mri_spim_reg:
    input:
        inputs["spim"].expand(
            bids(
                root=root,
                datatype="warps",
                space="SPIM",
                desc="deformwarped",
                suffix="T2w.nii.gz",
                iters="{iters}",
                dof="{dof}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["spim"].wildcards,
            ),
            iters="100x50x50x0",
            gradsigma=range(2, 4),
            warpsigma=range(2, 4),
            dof=[12],
            radius=[f"{i}x{i}x{i}" for i in range(2, 4)],
        ),
    group:
        "subj"


rule warp_mri_to_template_via_spim:
    input:
        mri=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards,
        ),
        ref=rules.import_template_anat.output.anat,
        affine_mri_to_spim=bids(
            root=root,
            datatype="warps",
            from_="mri",
            to="spim",
            type_="ras",
            desc="rigid",
            suffix="xfm.txt",
            iters="100x100x50x0",
            dof="12",
            radius="2x2x2",
            gradsigma="3",
            warpsigma="3",
            **inputs["spim"].wildcards,
        ),
        warp_mri_to_spim=bids(
            root=root,
            datatype="warps",
            from_="mri",
            to="spim",
            suffix="warp.nii.gz",
            iters="100x100x50x0",
            dof="12",
            radius="2x2x2",
            gradsigma="3",
            warpsigma="3",
            **inputs["spim"].wildcards,
        ),
        affine_spim_to_template=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
        warp_spim_to_template=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        warped=bids(
            root=root,
            datatype="anat",
            space="{template}",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        "  -rm {input.mri} {output.warped} "
        "  -r {input.warp_spim_to_template} {input.affine_spim_to_template} {input.warp_mri_to_spim} {input.affine_mri_to_spim}"


rule warp_mri_brainmask_to_spim:
    """Warp MRI brain mask to SPIM space for tissue comparison.
    
    Applies the MRI->SPIM transformation to the MRI brain mask to enable
    direct comparison of brain tissue extent between modalities. Also computes
    Jacobian determinant to quantify local tissue deformation due to fixation
    and clearing processes.
    """
    input:
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            iters="100x100x50",
            radius="2x2x2",
            gradsigma="3",
            warpsigma="3",
            **inputs["T2w"].wildcards,
        ),
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        affine_mri_to_spim=bids(
            root=root,
            datatype="warps",
            from_="mri",
            to="spim",
            type_="ras",
            desc="rigid",
            suffix="xfm.txt",
            iters="100x100x50x0",
            dof="12",
            radius="2x2x2",
            gradsigma="3",
            warpsigma="3",
            **inputs["spim"].wildcards,
        ),
        warp_mri_to_spim=bids(
            root=root,
            datatype="warps",
            from_="mri",
            to="spim",
            suffix="warp.nii.gz",
            iters="100x100x50x0",
            dof="12",
            radius="2x2x2",
            gradsigma="3",
            warpsigma="3",
            **inputs["spim"].wildcards,
        ),
    output:
        mask=bids(
            root=root,
            datatype="anat",
            space="spim",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
        jacobian=bids(
            root=root,
            datatype="anat",
            space="spim",
            desc="brain",
            suffix="jacobian.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=15,
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} -ri NN"
        "  -rm {input.mask} {output.mask} "
        "  -r {input.warp_mri_to_spim} {input.affine_mri_to_spim} -rj {output.jacobian}"
