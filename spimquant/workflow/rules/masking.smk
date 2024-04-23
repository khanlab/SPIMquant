rule atropos_seg:
    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    params:
        downsampling="50%",
        mrf_smoothing=0.3,
        mrf_radius="2x2x2",
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="Atropos",
            k="{k}",
            suffix="dseg.nii",
            **inputs["spim"].wildcards
        ),
        posteriors_dir=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="Atropos",
                k="{k}",
                suffix="posteriors",
                **inputs["spim"].wildcards
            )
        ),
    container:
        None
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=16000,
    shell:
        "mkdir -p {output.posteriors_dir} && "
        "c3d {input.nii} -resample {params.downsampling} -o downsampled.nii -scale 0 -shift 1 -o ones.nii && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "Atropos -v -d 3 --initialization KMeans[{wildcards.k}] "
        " --intensity-image downsampled.nii "
        " --output [dseg_downsampled.nii,{output.posteriors_dir}/class-%02d.nii] "
        " --mask-image ones.nii --mrf [{params.mrf_smoothing},{params.mrf_radius}] && "
        "c3d -interpolation NearestNeighbor {input.nii} dseg_downsampled.nii -reslice-identity -o {output.dseg}"


rule init_affine_reg:
    """initial affine registration used to obtain priors for brainmasking"""
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=config["masking"]["stain"],
            level=config["masking"]["level"],
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
            desc="initaffine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="{template}",
            desc="initaffinewarped",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
    log:
        bids(
            root="logs",
            datatype="init_affine_reg",
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


rule transform_template_dseg_to_subject:
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=config["masking"]["stain"],
            level=config["masking"]["level"],
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        dseg=rules.import_dseg.output.dseg,
        xfm_ras=rules.init_affine_reg.output.xfm_ras,
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            desc="initaffine",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards
        ),
    threads: 32
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        " -ri NN"
        "  -rm {input.dseg} {output.dseg} "
        "  -r {input.xfm_ras},-1 "


rule create_mask_from_gmm_and_prior:
    input:
        tissue_dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="Atropos",
            k=config["masking"]["gmm_k"],
            suffix="dseg.nii",
            **inputs["spim"].wildcards
        ),
        atlas_dseg=bids(
            root=root,
            datatype="micr",
            desc="initaffine",
            from_=config["masking"]["atlas_prior_for_mask"],
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards
        ),
    params:
        k=config["masking"]["gmm_k"],
    output:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),
    script:
        "../scripts/create_mask_from_gmm_and_prior.py"


rule create_mask_from_gmm:
    input:
        dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="Atropos",
            k=config["masking"]["gmm_k"],
            suffix="dseg.nii",
            **inputs["spim"].wildcards
        ),
    params:
        bg_label=config["masking"]["gmm_bg_class"],
    output:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain1class",
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),
    shell:
        "c3d {input} -threshold {params.bg_label} {params.bg_label} 0 1 -o {output}"
