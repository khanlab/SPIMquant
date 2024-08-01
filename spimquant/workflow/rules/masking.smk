rule pre_atropos:
    """Obtain down-sampled scanning image and a mask of ones of the same size
    
    This rule takes the down-sampled brain scan image as input, further down-sample it and also output a mask of 
    the same size
    
    input:
        nii: down-sampled brain scan image in .nii format (float32)
    
    output:
        downsampled: The further down-sampled brain scan image in .nii format (float32)
        mask: An image of ones of the same size as `downsampled` (float32)
    """

    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii",
            **cconfig.inputs["spim"].wildcards
        ),
    params:
        downsampling=config["masking"]["pre_atropos_downsampling"],
    output:
        downsampled=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="preAtropos",
                suffix="SPIM.nii",
                **cconfig.inputs["spim"].wildcards
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
                **cconfig.inputs["spim"].wildcards
            )
        ),
    container:
        config["containers"]["itksnap"]
    shell:
        "c3d {input.nii} -resample {params.downsampling} -o {output.downsampled} -scale 0 -shift 1 -o {output.mask}"


rule atropos_seg:
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
                **cconfig.inputs["spim"].wildcards
            )
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
                **cconfig.inputs["spim"].wildcards
            )
        ),
    container:
        config["containers"]["ants"]
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=16000,
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
            suffix="SPIM.nii",
            **cconfig.inputs["spim"].wildcards
        ),
        dseg=rules.atropos_seg.output.dseg,
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="Atropos",
            k="{k}",
            suffix="dseg.nii",
            **cconfig.inputs["spim"].wildcards
        ),
    container:
        config["containers"]["itksnap"]
    shell:
        "c3d -interpolation NearestNeighbor {input.ref} {input.dseg} -reslice-identity -o {output.dseg}"


rule init_affine_reg:
    """initial affine registration used to obtain priors for brainmasking
    
    This step take the observed brain scan data, and the template image and tries to find an affine transformation that 
    maps between the two (since they are of different size and scale).
    
    The affine transform obtained in this step is used to back transform templates onto observed subject brain scan.
    
    input:
        template: anat.nii image (uint8)
        subject: Down-sampled brain image (float64)
    
    output:
        xfm_ras: The affine transformation as .txt file
        warped: (float64) The warped (and downsampled) .nii brain scan onto the size of template image provided
    """
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=config["masking"]["stain"],
            level=config["masking"]["level"],
            suffix="SPIM.nii",
            **cconfig.inputs["spim"].wildcards
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
            **cconfig.inputs["spim"].wildcards
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="{template}",
            desc="initaffinewarped",
            suffix="SPIM.nii",
            **cconfig.inputs["spim"].wildcards
        ),
    log:
        bids(
            root="logs",
            datatype="init_affine_reg",
            space="{template}",
            suffix="log.txt",
            **cconfig.inputs["spim"].wildcards
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


rule affine_transform_template_dseg_to_subject:
    """
    The result of this rule is a downsampled .nii.gz image that is roughly matched against the OME zarr array. This 
    array is (as with dseg image) an instance segmentation mask of brain regions.
    """
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=config["masking"]["stain"],
            level=config["masking"]["level"],
            suffix="SPIM.nii",
            **cconfig.inputs["spim"].wildcards
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
            **cconfig.inputs["spim"].wildcards
        ),
    threads: 32
    container:
        config["containers"]["itksnap"]
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        " -ri NN"
        "  -rm {input.dseg} {output.dseg} "
        "  -r {input.xfm_ras},-1 "


rule create_mask_from_gmm_and_prior:
    """Create a down-sampled mask of brain scan
    
    input:
        tissue_dseg: 
        atlas_dseg: warped template image from affine_transform_template_dseg_to_subject (float32)
    output:
        mask: A down-sampled 0-1 binary .nii mask (float32), where 0 are background and 1 are brain region
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
            **cconfig.inputs["spim"].wildcards
        ),
        atlas_dseg=bids(
            root=root,
            datatype="micr",
            desc="initaffine",
            from_=config["masking"]["priors_template"],
            suffix="dseg.nii.gz",
            **cconfig.inputs["spim"].wildcards
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
            **cconfig.inputs["spim"].wildcards
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
            **cconfig.inputs["spim"].wildcards
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
            **cconfig.inputs["spim"].wildcards
        ),
    container:
        config["containers"]["itksnap"]
    shell:
        "c3d {input} -threshold {params.bg_label} {params.bg_label} 0 1 -o {output}"
