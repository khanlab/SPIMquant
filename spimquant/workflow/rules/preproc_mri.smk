rule import_subj_T2w:
    #    """ TODO: this sets orientation too, testing out L-R flip.."""
    input:
        nii=inputs["T2w"].path,
    output:
        nii=bids(
            root=root,
            datatype="anat",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards
        ),
    container: 
        config["containers"]["itksnap"]
    shell:
        #        'c3d {input} -orient RSP -o {output}' 
        "cp {input} {output}"



rule n4_mri:
    input:
        nii=bids(
            root=root, datatype="anat", suffix="T2w.nii.gz", **inputs["T2w"].wildcards
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards
        ),
    container:
        config["containers"]["ants"]
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
            **inputs["T2w"].wildcards
        ),
    params:
        iters="{iters}",#"100x50x50",
        metric="NCC {radius}",#3x3x3",
        sigma1="{gradsigma}vox", #3
        sigma2="{warpsigma}vox", #2
    output:
        xfm_ras=temp(bids(
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
            **inputs["T2w"].wildcards
        )),
        warp=temp(bids(
            root=root,
            datatype="warps",
            from_="mri",
            to="{template}",
            suffix="warp.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["T2w"].wildcards
        )),
        invwarp=temp(bids(
            root=root,
            datatype="warps",
            from_="{template}",
            to="mri",
            suffix="warp.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["T2w"].wildcards
        )),
        warped=temp(bids(
            root=root,
            datatype="warps",
            space="{template}",
            desc="deformwarped",
            suffix="T2w.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs["T2w"].wildcards
        )),
#    log:
#        bids(
#            root="logs",
#            datatype="mri_rigid_greedy",
#            space="{template}",
#            suffix="log.txt",
#            **inputs["T2w"].wildcards
#        ),
    threads: 32
    container:
        config["containers"]["itksnap"]
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
        expand(bids(
            root=root,
            datatype="anat",
            desc="brain",
            suffix="mask.nii.gz",
            iters="{iters}",
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",
            **inputs["T2w"].wildcards
        ),iters="100x50x50",
               gradsigma=range(3,6),
               warpsigma=range(3,6),
               radius=[f"{i}x{i}x{i}" for i in range(2,5)],
               subject='o21')





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
            **inputs["T2w"].wildcards
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
            **inputs["T2w"].wildcards
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
            **inputs["T2w"].wildcards
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
            **inputs["T2w"].wildcards
        ),
    shadow:
        "minimal"
    threads: 32
    container:
        config["containers"]["itksnap"]
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
            **inputs["T2w"].wildcards
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
            **inputs["T2w"].wildcards
        ),
    output:
        nii=bids(
            root=root,
            datatype="anat",
            desc="N4brain",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards
        ),
    container:
        config["containers"]["itksnap"]
    shell:
        "c3d {input.nii} {input.mask} -multiply -resample 300% -o {output.nii}"


rule rigid_greedy_reg_mri_to_spim:
    input:
        mri=bids(
            root=root,
            datatype="anat",
            desc="N4brain",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards
        ),
        spim=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    params:
        iters="{iters}",#"100x50x50",
        metric_rigid="NMI",
        dof='{dof}',
        metric="NCC {radius}",#3x3x3",
        sigma1="{gradsigma}vox", #3
        sigma2="{warpsigma}vox", #2
 
    output:
        xfm_ras=temp(bids(
            root=root,
            datatype="warps",
            from_="mri",
            to="spim",
            type_="ras",
            desc="rigid",
            suffix="xfm.txt",
            iters="{iters}",
            dof='{dof}',
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["spim"].wildcards
        )),
        warp=temp(bids(
            root=root,
            datatype="warps",
            from_="mri",
            to="spim",
            suffix="warp.nii.gz",
            iters="{iters}",
            dof='{dof}',
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["spim"].wildcards
        )),
        invwarp=temp(bids(
            root=root,
            datatype="warps",
            from_="spim",
            to="mri",
            suffix="warp.nii.gz",
            iters="{iters}",
            dof='{dof}',
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["spim"].wildcards
        )),
        linwarped=bids(
            root=root,
            datatype="warps",
            space="SPIM",
            desc="linearwarped",
            suffix="T2w.nii.gz",
            iters="{iters}",
            dof='{dof}',
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["spim"].wildcards
        ),

        warped=bids(
            root=root,
            datatype="warps",
            space="SPIM",
            desc="deformwarped",
            suffix="T2w.nii.gz",
            iters="{iters}",
            dof='{dof}',
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["spim"].wildcards
        ),
#    log:
#        bids(
#            root="logs",
#            datatype="mri_spim_rigid_greedy",
#            suffix="log.txt",
#            **inputs["spim"].wildcards
#        ),
    threads: 32
    container:
        config["containers"]["itksnap"]
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
        expand(bids(
            root=root,
            datatype="warps",
            space="SPIM",
            desc="deformwarped",
            suffix="T2w.nii.gz",
            iters="{iters}",
            dof='{dof}',
            radius="{radius}",
            gradsigma="{gradsigma}",
            warpsigma="{warpsigma}",

            **inputs["spim"].wildcards
        ),iters="100x50x50x0",
               gradsigma=range(2,4),
               warpsigma=range(2,4),
               dof=[12],
               radius=[f"{i}x{i}x{i}" for i in range(2,4)],
               subject='o21',
               sample='brain',
               acq='blaze')

rule warp_mri_to_template_via_spim:
    input:
        mri=bids(
            root=root,
            datatype="anat",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["T2w"].wildcards
        ),
        ref=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"),
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
            **inputs["spim"].wildcards
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
            **inputs["spim"].wildcards
        ),
        affine_spim_to_template=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards
        ),
        warp_spim_to_template=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii",
            **inputs["spim"].wildcards
        ),

    output:
        warped=bids(
            root=root,
            datatype="anat",
            space="{template}",
            desc="N4",
            suffix="T2w.nii.gz",
            **inputs["spim"].wildcards
        ),

    threads: 32
    container:
        config["containers"]["itksnap"]
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
            **inputs["T2w"].wildcards
        ),
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            desc=config["templatereg"]["desc"],
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
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
            **inputs["spim"].wildcards
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
            **inputs["spim"].wildcards
        ),
    output:
        mask=bids(
            root=root,
            datatype="anat",
            space="spim",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards
        ),
        jacobian=bids(
            root=root,
            datatype="anat",
            space="spim",
            desc="brain",
            suffix="jacobian.nii.gz",
            **inputs["spim"].wildcards
        ),

    threads: 32
    container:
        config["containers"]["itksnap"]
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} -ri NN"
        "  -rm {input.mask} {output.mask} "
        "  -r {input.warp_mri_to_spim} {input.affine_mri_to_spim} -rj {output.jacobian}"




