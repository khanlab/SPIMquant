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


rule n4_mri:
    input:
        nii=select_single_mri,
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


rule greedy_reg_mri_to_spim:
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
    """ to assess effect of perfusion fixation and clearing"""
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
