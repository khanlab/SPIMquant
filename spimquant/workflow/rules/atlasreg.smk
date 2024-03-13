wildcard_constraints:
    level='[0-9]+'

rule get_downsampled_nii:
    input:
        zarr=inputs["spim"].path,
    params:
        channel_index=lambda wildcards: config["channel_mapping"][wildcards.stain],
    output:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="spim.nii",
            **inputs["spim"].wildcards
        ),
    threads: 32
    script:
        "../scripts/ome_zarr_to_nii.py"


rule import_anat:
    input:
        anat=lambda wildcards: format(config["atlases"][wildcards.template]["anat"]),
    output:
        anat=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    log:
        bids_tpl(
            root="logs",
            datatype="import_anat",
            template="{template}",
            suffix="log.txt",
        ),
    shell:
        "cp {input} {output}"


rule import_dseg:
    input:
        dseg=lambda wildcards: format(config["atlases"][wildcards.template]["dseg"]),
    output:
        dseg=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
    log:
        bids_tpl(
            root="logs",
            datatype="import_dseg",
            template="{template}",
            suffix="log.txt",
        ),
    shell:
        "cp {input} {output}"


rule import_labelmapper_lut:
    input:
        json=lambda wildcards: format(config["atlases"][wildcards.template]["lut"]),
    output:
        tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
    log:
        bids_tpl(
            root="logs", datatype="import_lut", template="{template}", suffix="log.txt"
        ),
    script:
        "../scripts/import_labelmapper_lut.py"


    


rule make_itksnap_lut:
    input:
        tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
    output:
        lut=bids_tpl(root=root, template="{template}", desc='itksnap',suffix="labels.txt"),
    script:
        "../scripts/lut_bids_to_itksnap.py"



# --- processing to get mask
rule atropos_seg:
    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="spim.nii",
            **inputs["spim"].wildcards
        ),
    params:
        mrf_smoothing=0.3,
        mrf_radius='2x2x2',
    output:
        dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='Atropos',
            k='{k}',
            suffix="dseg.nii",
            **inputs["spim"].wildcards
        ),
        posteriors_dir=directory(bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='Atropos',
            k='{k}',
            suffix="posteriors",
            **inputs["spim"].wildcards
        )),

    container:
        None
    shadow: 'minimal'
    shell:
        "mkdir -p {output.posteriors_dir} && "
        "c3d {input.nii} -scale 0 -shift 1 ones.nii && "
        "Atropos -v -d 3 --initialization KMeans[{wildcards.k}] "
        " --intensity-image {input.nii} "
        " --output [{output.dseg},{output.posteriors_dir}/class-%02d.nii] "
        " --mask-image ones.nii --mrf [{params.mrf_smoothing},{params.mrf_radius}]"

rule init_affine_reg:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=config["atlasreg"]["stain"],
            level=config["atlasreg"]["level"],
            suffix="spim.nii",
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
            suffix="spim.nii",
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
    shell:
        "greedy -d 3 -i {input.template} {input.subject} "
        " -a -dof 12 -ia-image-centers -m NMI -o {output.xfm_ras} && "
        " greedy -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.xfm_ras}"

rule transform_template_dseg_to_subject:
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=config["atlasreg"]["stain"],
            level=config["atlasreg"]["level"],
            suffix="spim.nii",
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
            )
    shell:
        " greedy -d 3 -rf {input.ref} "
        "  -rm {input.dseg} {output.dseg} "
        "  -r {input.xfm_ras},-1 "
        " -ri NN"

rule create_mask_from_gmm_and_prior:
    input:
        tissue_dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='Atropos',
            k=config['masking']['gmm_k'],
            suffix="dseg.nii",
            **inputs["spim"].wildcards
        ),
        atlas_dseg=bids(
                    root=root,
                    datatype="micr",
                    desc="initaffine",
                    from_=config['atlasreg']['atlas_prior_for_mask'],
                    suffix="dseg.nii.gz",
                    **inputs["spim"].wildcards
            )
    params:
            k=config['masking']['gmm_k'],
    output:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='brain',
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),  
    script:
        '../scripts/create_mask_from_gmm_and_prior.py'
    

rule create_mask_from_gmm:
    input:
        dseg=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='Atropos',
            k=config['masking']['gmm_k'],
            suffix="dseg.nii",
            **inputs["spim"].wildcards
        ),
    params:
        bg_label=config['masking']['gmm_bg_class']
    output:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='brain1class',
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),  
    shell:
        'c3d {input} -threshold {params.bg_label} {params.bg_label} 0 1 -o {output}'

    
rule n4:
    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="spim.nii",
            **inputs["spim"].wildcards
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='brain',
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),  

    output:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='N4',
            suffix="spim.nii",
            **inputs["spim"].wildcards
        ),
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='N4',
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
            desc='N4',
            suffix="spim.nii",
            **inputs["spim"].wildcards
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='brain',
            suffix="mask.nii",
            **inputs["spim"].wildcards)
    output:
        masked=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc='N4brain',
            suffix="spim.nii",
            **inputs["spim"].wildcards
        ),
    shell:
        'c3d {input.corrected} {input.mask} -multiply -o {output.masked}'

rule affine_reg:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            stain=config["atlasreg"]["stain"],
            level=config["atlasreg"]["level"],
            desc=config["atlasreg"]["desc"],
            suffix="spim.nii",
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
            suffix="spim.nii",
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
    shell:
        "greedy -d 3 -i {input.template} {input.subject} "
        " -a -dof 12 -ia-image-centers -m NMI -o {output.xfm_ras} && "
        " greedy -d 3 -rf {input.template} "
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
            suffix="spim.nii",
            **inputs["spim"].wildcards
        ),
        xfm_ras=rules.affine_reg.output.xfm_ras,
    params:
        iters='100x50',
        metric='NMI',
        sigma1='4vox',
        sigma2='2vox'
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
            suffix="spim.nii",
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
    shell:
        "greedy -d 3 -i {input.template} {input.subject} "
        " -it {input.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} && "
        " greedy -d 3 -rf {input.template} "
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

rule affine_to_template_nii:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        channel_index=lambda wildcards: config["channel_mapping"][wildcards.stain],
    output:
        nii=bids(
                    root=root,
                    datatype="micr",
                    desc="affine",
                    space="{template}",
                    stain="{stain}",
                    suffix="spim.nii",
                    **inputs["spim"].wildcards
                )
    container: None
    threads: 32
    script: '../scripts/affine_to_template_nii.py'

rule affine_to_template_ome_zarr:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        channel_index=lambda wildcards: config["channel_mapping"][wildcards.stain],
    output:
        ome_zarr=directory(bids(
                    root=root,
                    datatype="micr",
                    desc="affine",
                    space="{template}",
                    stain="{stain}",
                    suffix="spim.ome.zarr",
                    **inputs["spim"].wildcards
                ))
    container: None
    threads: 32
    script: '../scripts/affine_to_template_ome_zarr.py'

rule deform_to_template_nii:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        channel_index=lambda wildcards: config["channel_mapping"][wildcards.stain],
        chunks=(50,50,50),
        zooms=None,
    output:
        nii=bids(
                    root=root,
                    datatype="micr",
                    desc="deform",
                    space="{template}",
                    stain="{stain}",
                    suffix="spim.nii",
                    **inputs["spim"].wildcards
                )
    container: None
    threads: 32
    script: '../scripts/deform_to_template_nii.py'


rule deform_to_template_nii_zoomed:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        channel_index=lambda wildcards: config["channel_mapping"][wildcards.stain],
        chunks=(50,50,50),
        zooms=lambda wildcards: (float(wildcards.res)/1000,float(wildcards.res)/1000,float(wildcards.res)/1000) #None #same resolution as template if NOne
    output:
        nii=bids(
                    root=root,
                    datatype="micr",
                    desc="deform",
                    space="{template}",
                    stain="{stain}",
                    res="{res}um",
                    suffix="spim.nii",
                    **inputs["spim"].wildcards
                )
    container: None
    threads: 32
    script: '../scripts/deform_to_template_nii.py'

rule deform_to_template_nii_nb:
    input:
        ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        warp_nii=rules.deform_reg.output.warp,
        ref_nii=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
    params:
        channel_index=lambda wildcards: config["channel_mapping"][wildcards.stain],
        chunks=(20,20,20),
        zooms=None #same resolution as template if NOne
    output:
        nii=bids(
                    root=root,
                    datatype="micr",
                    desc="deformnb",
                    space="{template}",
                    stain="{stain}",
                    suffix="spim.nii",
                    **inputs["spim"].wildcards
                )
    container: None
    threads: 32
    notebook: '../notebooks/deform_to_template_nii.py.ipynb'


rule deform_transform_labels_to_subj:
    input:
        ref_ome_zarr=inputs["spim"].path,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        invwarp_nii=rules.deform_reg.output.invwarp,
        flo_nii=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
    output:
        zarr=directory(bids(
                    root=root,
                    datatype="micr",
                    desc="deform",
                    space="subject",
                    suffix="dseg.zarr",
                    **inputs["spim"].wildcards
                ))
    container: None
    threads: 32
    script: '../scripts/deform_transform_channel_to_template_nii.py'



rule transform_labels_to_zoomed_template:
    input:
        dseg=bids_tpl(root=root, template="{template}", suffix="dseg.nii.gz"),
        ref=bids(
                    root=root,
                    datatype="micr",
                    desc="deform",
                    space="{template}",
                    stain=config['atlasreg']['stain'],
                    res="{res}um",
                    suffix="spim.nii",
                    **inputs["spim"].wildcards
                )
    output:
        dseg=bids(
                    root=root,
                    datatype="micr",
                    space="{template}",
                    res="{res}um",
                    suffix="dseg.nii",
                    **inputs["spim"].wildcards
                )
    shell:
        "antsApplyTransforms -d 3 -v -n NearestNeighbor "
        " -i {input.dseg} -o {output.dseg} "
        " -r {input.ref} "

rule ome_zarr_to_zipstore:
    """ generic rule to process any ome.zarr from work """
    input:
        zarr=f"{work}/{{prefix}}.ome.zarr",
    output:
        zarr_zip=f"{root}/{{prefix}}.ome.zarr.zip",
    log:
        "logs/ome_zarr_to_zipstore/{prefix}.log",
    group:
        "preproc"
    shell:
        "7z a -mx0 -tzip {output.zarr_zip} {input.zarr}/. &> {log}"


