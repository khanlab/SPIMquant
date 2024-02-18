#adjust these to use snakebids: inputs['spim'].wildcards, 
print('inputs:')       
print(inputs['spim'])
print('outputs:')
print(bids(root=root,datatype='micr',desc='downsampled',suffix='spim.nii',**inputs['spim'].wildcards))


rule get_downsampled_nii:
    input:
        zarr=inputs['spim'].path
    params:
        level=config['atlasreg']['level']
    output:
        nii=bids(
            root=root,
            datatype="micr",
            suffix="spim.nii",
            **inputs['spim'].wildcards
        ),
    threads: 32
    script:
        "../scripts/ome_zarr_to_nii.py"

rule import_anat:
    input:
        anat=lambda wildcards: config["atlases"][wildcards.template]["anat"],
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
        dseg=lambda wildcards: config["atlases"][wildcards.template]["dseg"],
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


rule import_lut:
    input:
        json=lambda wildcards: config["atlases"][wildcards.template]["lut"],
    output:
        tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
    log:
        bids_tpl(
            root="logs", datatype="import_lut", template="{template}", suffix="log.txt"
        ),
    script:
        "../scripts/import_labelmapper_lut.py"


rule affine_reg:
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        subject=bids(
            root=root,
            datatype="micr",
            suffix="spim.nii",
            **inputs['spim'].wildcards
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
            **inputs['spim'].wildcards
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="{template}",
            desc="affinewarped",
            suffix="spim.nii",
            **inputs['spim'].wildcards
        ),
    log:
        bids(
            root="logs",
            datatype="affine_reg",
            space="{template}",
            suffix="log.txt",
            **inputs['spim'].wildcards
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
            suffix="spim.nii",
            **inputs['spim'].wildcards
        ),
        xfm_ras=rules.affine_reg.output.xfm_ras,
    output:
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii",
            **inputs['spim'].wildcards
        ),
        warped=bids(
            root=root,
            datatype="warps",
            space="{template}",
            desc="deformwarped",
            suffix="spim.nii",
            **inputs['spim'].wildcards
        ),
    log:
        bids(
            root="logs",
            datatype="deform_reg",
            space="{template}",
            suffix="log.txt",
            **inputs['spim'].wildcards
        ),
    shell:
        "greedy -d 3 -i {input.template} {input.subject} "
        " -it {input.xfm_ras} -m NMI "
        " -o {output.warp} -n 100x50x0x0 && "
        " greedy -d 3 -rf {input.template} "
        "  -rm {input.subject} {output.warped} "
        "  -r {output.warp} {input.xfm_ras}"


rule resample_labels_to_zarr:
    """TODO: add required OME metadata"""
    input:
        dseg=rules.import_dseg.output.dseg,
        xfm_ras=rules.affine_reg.output.xfm_ras,
        label_tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
        zarr_zip=inputs['spim'].path
    params:
        level_to_resample_to=0,
        max_downsampling_layers=config['ome_zarr']['max_downsampling_layers'],
        label_name='dseg',
        scaling_method='nearest'
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    desc="resampled",
                    from_="{template}",
                    suffix="dseg.ome.zarr",
                    **inputs['spim'].wildcards
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
            **inputs['spim'].wildcards
        ),
    script:
        "../scripts/resample_labels_to_zarr.py"


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


