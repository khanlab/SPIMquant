"""
MRI preprocessing and cross-modality registration workflow for SPIMquant.

This module handles co-registration of in-vivo MRI (T2w) with ex-vivo SPIM data,
enabling multi-modal analysis and assessment of tissue changes due to perfusion
fixation and optical clearing.

Key workflow stages:
1. N4 bias field correction of MRI
2. MRI to MRI-template rigid+nlin registration (for brain masking)
3. Template brain mask to MRI transformation
4. MRI brain extraction
5. MRI to SPIM affine+nlin registration
6. Parameter tuning rules for optimization
7. Concatenated transformations (MRI -> SPIM -> Template)
8. Jacobian determinant calculation (tissue deformation quantification)
9. Registration QC report generation

This workflow is optional and used when both MRI and SPIM data are available
for the same subject. 
"""


def select_single_mri(wildcards):

    files = inputs["mri"].filter(subject=wildcards.subject).expand()
    if len(files) == 1:
        return files[0]
    elif len(files) == 0:
        raise ValueError(f"No MRI images found for f{wildcards}")
    else:
        raise ValueError(
            f"Multiple MRI images found for f{wildcards}, use --filter-mri to select a single image"
        )


def get_all_mri(wildcards):
    """Get all MRI files for a subject (for multi-MRI averaging)."""
    return inputs["mri"].filter(subject=wildcards.subject).expand()


def get_mri_indices(wildcards):
    """Get list of indices for MRI files for a subject."""
    files = inputs["mri"].filter(subject=wildcards.subject).expand()
    return list(range(len(files)))


rule n4_mri_individual:
    """Apply N4 bias field correction to individual T2w MRI images.
    
    Uses ANTs N4BiasFieldCorrection to correct intensity inhomogeneities in
    each MRI image before registration and averaging.
    """
    input:
        nii=lambda wildcards: get_all_mri(wildcards)[int(wildcards.mriidx)],
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4",
            mriidx="{mriidx}",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
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


rule register_mri_to_first:
    """Rigidly register each MRI to the first MRI using greedy.
    
    This aligns multiple MRI acquisitions to enable super-resolved averaging.
    Only runs when multiple MRI images are present.
    """
    input:
        fixed=bids(
            root=root,
            datatype="anat",
            desc="N4",
            mriidx="0",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        moving=bids(
            root=root,
            datatype="anat",
            desc="N4",
            mriidx="{mriidx}",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    output:
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_=f"{mri_suffix}",
            to=f"{mri_suffix}ref",
            type_="ras",
            desc="rigid",
            mriidx="{mriidx}",
            suffix="xfm.txt",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=10,
    shell:
        "if [ {wildcards.mriidx} -eq 0 ]; then "
        "  echo '1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1' > {output.xfm_ras}; "
        "else "
        "  greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} "
        "  -a -dof 6 -ia-image-centers -m NMI -o {output.xfm_ras}; "
        "fi"


rule average_mri:
    """Average all registered and N4-corrected MRI images.
    
    Optionally upsamples images during transformation for super-resolution.
    Produces a final preprocessed MRI with desc-preproc.
    """
    input:
        ref=bids(
            root=root,
            datatype="anat",
            desc="N4",
            mriidx="0",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        mri=lambda wildcards: expand(
            bids(
                root=root,
                datatype="anat",
                desc="N4",
                mriidx="{mriidx}",
                suffix=f"{mri_suffix}.nii.gz",
                **inputs.subj_wildcards,
            ),
            mriidx=get_mri_indices(wildcards),
            allow_missing=True,
        ),
        xfm=lambda wildcards: expand(
            bids(
                root=root,
                datatype="warps",
                from_=f"{mri_suffix}",
                to=f"{mri_suffix}ref",
                type_="ras",
                desc="rigid",
                mriidx="{mriidx}",
                suffix="xfm.txt",
                **inputs.subj_wildcards,
            ),
            mriidx=get_mri_indices(wildcards),
            allow_missing=True,
        ),
    params:
        upsample_factor=lambda wildcards: config.get("mri_upsample_factor", 1),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="preproc",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/average_registered_mri.py"


rule n4_mri:
    """Apply N4 bias field correction to T2w MRI.
    
    Uses ANTs N4BiasFieldCorrection to correct intensity inhomogeneities in
    the MRI image, improving subsequent registration performance.
    
    This rule handles both single and multiple MRI cases. For multiple MRIs,
    it uses the preprocessed (averaged) output which is already N4-corrected.
    For single MRIs, it applies N4 directly.
    """
    input:
        nii=lambda wildcards: (
            bids(
                root=root,
                datatype="anat",
                desc="preproc",
                suffix=f"{mri_suffix}.nii.gz",
                **inputs.subj_wildcards,
            )
            if len(get_all_mri(wildcards)) > 1
            else select_single_mri(wildcards)
        ),
    params:
        is_multi_mri=lambda wildcards: len(get_all_mri(wildcards)) > 1,
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    script:
        "../scripts/copy_or_n4.py"


rule rigid_nlin_reg_mri_to_template:
    """Initial unmasked MRI to unmasked template MRI using rigid + 
    deformable registration.
    
    Performs initial rigid (6 DOF) alignment followed by deformable
    registration to align non-brain-masked  MRI with similarly
    non-brain-masked MRI, so that the template brain mask can be
    propagated. Note: this step could be replaced by a ML-based 
    brain-masking if a suitable model is found/trained.
    """
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
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
                from_=f"{mri_suffix}",
                to="{template}",
                type_="ras",
                desc="rigid",
                suffix="xfm.txt",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
        warp=temp(
            bids(
                root=root,
                datatype="warps",
                from_=f"{mri_suffix}",
                to="{template}",
                suffix="warp.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
        invwarp=temp(
            bids(
                root=root,
                datatype="warps",
                from_="{template}",
                to=f"{mri_suffix}",
                suffix="warp.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
            )
        ),
        warped=temp(
            bids(
                root=root,
                datatype="warps",
                space="{template}",
                desc="deformwarped",
                suffix=f"{mri_suffix}.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
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


rule all_tune_mri_mask:
    input:
        inputs["mri"].expand(
            bids(
                root=root,
                datatype="anat",
                desc="brain",
                suffix="mask.nii.gz",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs.subj_wildcards,
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
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_=f"{mri_suffix}",
            to=config["template_mri"],
            type_="ras",
            desc="rigid",
            suffix="xfm.txt",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs.subj_wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            from_=config["template_mri"],
            to=f"{mri_suffix}",
            suffix="warp.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs.subj_wildcards,
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
            **inputs.subj_wildcards,
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
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        mask=bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            iters=config["reg_mri"]["rigidgreedy"]["iters"],
            radius=config["reg_mri"]["rigidgreedy"]["radius"],
            gradsigma=config["reg_mri"]["rigidgreedy"]["gradsigma"],
            warpsigma=config["reg_mri"]["rigidgreedy"]["warpsigma"],
            **inputs.subj_wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4brain",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
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


rule affine_nlin_reg_mri_to_spim:
    """Register MRI to SPIM space using affine + deformable registration.
    
    Performs initial affine (6 or 12 DOF) alignment followed by deformable
    registration to align brain-masked MRI with SPIM data. 
    """
    input:
        mri=bids(
            root=root,
            datatype="anat",
            desc="N4brain",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
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
                from_=f"{mri_suffix}",
                to="SPIM",
                type_="ras",
                desc="rigid",
                suffix="xfm.txt",
                dof="{dof}",
                iters="{iters}",
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
                from_=f"{mri_suffix}",
                to="SPIM",
                suffix="warp.nii.gz",
                dof="{dof}",
                iters="{iters}",
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
                from_="SPIM",
                to=f"{mri_suffix}",
                suffix="warp.nii.gz",
                dof="{dof}",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["spim"].wildcards,
            )
        ),
        linwarped=temp(
            bids(
                root=root,
                datatype="warps",
                space="SPIM",
                desc="linearwarped",
                suffix=f"{mri_suffix}.nii.gz",
                dof="{dof}",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["spim"].wildcards,
            )
        ),
        warped=temp(
            bids(
                root=root,
                datatype="warps",
                space="SPIM",
                desc="deformwarped",
                suffix=f"{mri_suffix}.nii.gz",
                dof="{dof}",
                iters="{iters}",
                radius="{radius}",
                gradsigma="{gradsigma}",
                warpsigma="{warpsigma}",
                **inputs["spim"].wildcards,
            )
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


rule all_tune_mri_spim_reg:
    input:
        inputs["spim"].expand(
            bids(
                root=root,
                datatype="warps",
                space="SPIM",
                desc="deformwarped",
                suffix=f"{mri_suffix}.nii.gz",
                dof="{dof}",
                iters="{iters}",
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
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        ref=rules.import_template_anat.output.anat,
        affine_mri_to_spim=bids(
            root=root,
            datatype="warps",
            from_=f"{mri_suffix}",
            to="SPIM",
            type_="ras",
            desc="rigid",
            suffix="xfm.txt",
            dof=config["reg_mri"]["greedy"]["dof"],
            iters=config["reg_mri"]["greedy"]["iters"],
            radius=config["reg_mri"]["greedy"]["radius"],
            gradsigma=config["reg_mri"]["greedy"]["gradsigma"],
            warpsigma=config["reg_mri"]["greedy"]["warpsigma"],
            **inputs["spim"].wildcards,
        ),
        warp_mri_to_spim=bids(
            root=root,
            datatype="warps",
            from_=f"{mri_suffix}",
            to="SPIM",
            suffix="warp.nii.gz",
            dof=config["reg_mri"]["greedy"]["dof"],
            iters=config["reg_mri"]["greedy"]["iters"],
            radius=config["reg_mri"]["greedy"]["radius"],
            gradsigma=config["reg_mri"]["greedy"]["gradsigma"],
            warpsigma=config["reg_mri"]["greedy"]["warpsigma"],
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
            via="SPIM",
            desc="N4",
            suffix=f"{mri_suffix}.nii.gz",
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
            **inputs.subj_wildcards,
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
            from_=f"{mri_suffix}",
            to="SPIM",
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
            from_=f"{mri_suffix}",
            to="SPIM",
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
            space="SPIM",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
        jacobian=bids(
            root=root,
            datatype="anat",
            space="SPIM",
            desc="brain",
            suffix="jacobian.nii.gz",
            **inputs["spim"].wildcards,
        ),
        composed_warp=bids(
            root=root,
            datatype="anat",
            from_=f"{mri_suffix}",
            to="SPIM",
            suffix="warp.nii.gz",
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


rule mri_spim_registration_qc_report:
    """Generate MRI to SPIM registration quality control report with visualizations"""
    input:
        spim=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        mri=bids(
            root=root,
            datatype="anat",
            desc="N4brain",
            suffix=f"{mri_suffix}.nii.gz",
            **inputs.subj_wildcards,
        ),
        warped_affine=bids(
            root=root,
            datatype="warps",
            space="SPIM",
            desc="linearwarped",
            suffix=f"{mri_suffix}.nii.gz",
            dof=config["reg_mri"]["greedy"]["dof"],
            iters=config["reg_mri"]["greedy"]["iters"],
            radius=config["reg_mri"]["greedy"]["radius"],
            gradsigma=config["reg_mri"]["greedy"]["gradsigma"],
            warpsigma=config["reg_mri"]["greedy"]["warpsigma"],
            **inputs["spim"].wildcards,
        ),
        warped_deform=bids(
            root=root,
            datatype="warps",
            space="SPIM",
            desc="deformwarped",
            suffix=f"{mri_suffix}.nii.gz",
            dof=config["reg_mri"]["greedy"]["dof"],
            iters=config["reg_mri"]["greedy"]["iters"],
            radius=config["reg_mri"]["greedy"]["radius"],
            gradsigma=config["reg_mri"]["greedy"]["gradsigma"],
            warpsigma=config["reg_mri"]["greedy"]["warpsigma"],
            **inputs["spim"].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            from_=f"{mri_suffix}",
            to="SPIM",
            suffix="warp.nii.gz",
            dof=config["reg_mri"]["greedy"]["dof"],
            iters=config["reg_mri"]["greedy"]["iters"],
            radius=config["reg_mri"]["greedy"]["radius"],
            gradsigma=config["reg_mri"]["greedy"]["gradsigma"],
            warpsigma=config["reg_mri"]["greedy"]["warpsigma"],
            **inputs["spim"].wildcards,
        ),
    params:
        stain_for_reg=stain_for_reg,
    output:
        report_html=bids(
            root=root,
            datatype="anat",
            space="SPIM",
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
        "../scripts/mri_spim_reg_qc_report.py"
