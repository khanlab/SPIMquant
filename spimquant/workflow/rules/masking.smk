"""
Brain masking workflow for SPIMquant.

This module creates brain tissue masks using Gaussian mixture modeling (Atropos)
combined with template-derived priors. The masks separate brain tissue from 
background and are used in subsequent intensity correction and registration steps.

Key workflow stages:
1. Downsampling and preprocessing for GMM segmentation
2. Atropos tissue classification (k-class GMM with MRF smoothing)
3. Initial affine registration to template
4. Transform template brain mask to subject space as prior
5. Combine GMM classes with template prior to create final brain mask

The workflow uses iterative registration to align a template-based prior with
tissue classification results, improving mask quality especially at brain boundaries.
"""


rule pre_atropos:
    """Prepare image for Atropos segmentation.
    
    Downsamples and preprocesses the SPIM image for efficient tissue classification.
    Applies log transformation and intensity normalization to improve GMM convergence.
    Also creates an initial foreground mask to restrict computation to brain regions.
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
    params:
        downsampling=(
            "10%"
            if config["sloppy"]
            else config["masking"]["pre_atropos_downsampling"]
        ),
    output:
        downsampled=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="preAtropos",
                suffix="SPIM.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
        mask=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="preAtropos",
                suffix="mask.nii",
                **inputs["spim"].wildcards,
            )
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
        "c3d {input.nii} -resample {params.downsampling} -shift 1 -log -stretch 2% 98% 0 100 -clip 0 100 -o {output.downsampled} -scale 0 -shift 1 -o {output.mask}"


rule atropos_seg:
    """Perform tissue classification using Atropos (k-class GMM).
    
    Uses ANTs Atropos to classify tissue into k intensity classes via Gaussian
    mixture modeling with Markov random field (MRF) spatial smoothing. Outputs
    a discrete segmentation and posterior probability maps for each class.
    """
    input:
        downsampled=rules.pre_atropos.output.downsampled,
        mask=rules.pre_atropos.output.mask,
    params:
        mrf_smoothing=0.3,
        mrf_radius="2x2x2",
    output:
        dseg=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="dsAtropos",
                k="{k}",
                suffix="dseg.nii",
                **inputs["spim"].wildcards,
            )
        ),
        posteriors_dir=temp(
            directory(
                bids(
                    root=root,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="Atropos",
                    k="{k}",
                    suffix="posteriors",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    conda:
        "../envs/ants.yaml"
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    shell:
        "mkdir -p {output.posteriors_dir} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "Atropos -v -d 3 --initialization KMeans[{wildcards.k}] "
        " --intensity-image {input.downsampled} "
        " --output [{output.dseg},{output.posteriors_dir}/class-%02d.nii] "
        " --mask-image {input.mask} --mrf [{params.mrf_smoothing},{params.mrf_radius}]"


rule post_atropos:
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        dseg=rules.atropos_seg.output.dseg,
    output:
        dseg=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="Atropos",
                k="{k}",
                suffix="dseg.nii",
                **inputs["spim"].wildcards,
            )
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
        "c3d -interpolation NearestNeighbor {input.ref} {input.dseg} -reslice-identity -o {output.dseg}"


rule init_affine_reg:
    """Perform initial affine registration for obtaining mask priors.
    
    Registers subject SPIM to template using 12-DOF affine transformation.
    This initial alignment enables template brain masks to be warped to subject
    space as priors for brain masking, even before final registration.
    """
    input:
        template=get_template_for_reg,
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    params:
        iters="10x0x0" if config["sloppy"] else "100x100",
    output:
        xfm_ras=temp(
            bids(
                root=root,
                datatype="warps",
                from_="subject",
                to="{template}",
                type_="ras",
                desc="initaffine",
                suffix="xfm.txt",
                **inputs["spim"].wildcards,
            )
        ),
        warped=temp(
            bids(
                root=root,
                datatype="warps",
                space="{template}",
                desc="initaffinewarped",
                suffix="SPIM.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj"
    log:
        bids(
            root="logs",
            datatype="init_affine_reg",
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


rule affine_transform_template_mask_to_subject:
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        dseg=rules.import_mask.output.mask,
        xfm_ras=rules.init_affine_reg.output.xfm_ras,
    output:
        dseg=temp(
            bids(
                root=root,
                datatype="micr",
                desc="initaffine",
                from_="{template}",
                suffix="mask.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        " -ri NN"
        "  -rm {input.dseg} {output.dseg} "
        "  -r {input.xfm_ras},-1 "


rule create_mask_from_gmm_and_prior:
    """Create final brain mask by combining GMM classes with template prior.
    
    Combines tissue classification results from Atropos with the template-derived
    brain mask to create a refined brain mask. Uses both intensity-based tissue
    classification and spatial prior information for improved accuracy.
    """
    input:
        tissue_dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="Atropos",
            k=config["masking"]["gmm_k"],
            suffix="dseg.nii",
            **inputs["spim"].wildcards,
        ),
        template_mask=bids(
            root=root,
            datatype="micr",
            desc="initaffine",
            from_=config["template"],
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
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
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
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
            **inputs["spim"].wildcards,
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
            suffix="mask.nii.gz",
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
        "c3d {input} -threshold {params.bg_label} {params.bg_label} 0 1 -o {output}"
