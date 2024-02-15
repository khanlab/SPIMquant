
#adjust these to use snakebids: inputs['spim'].wildcards, 

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
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc=config["atlasreg"]["desc"],
            stain=config["atlasreg"]["stain"],
            level=config["atlasreg"]["level"],
            suffix="spim.nii",
        ),
    output:
        xfm_ras=bids(
            root=root,
            subject="{subject}",
            datatype="warps",
            sample="{sample}",
            acq="{acq}",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
        ),
        warped=bids(
            root=root,
            subject="{subject}",
            datatype="warps",
            sample="{sample}",
            acq="{acq}",
            space="{template}",
            desc="affinewarped",
            suffix="spim.nii",
        ),
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="affine_reg",
            sample="{sample}",
            acq="{acq}",
            space="{template}",
            suffix="log.txt",
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
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc=config["atlasreg"]["desc"],
            stain=config["atlasreg"]["stain"],
            level=config["atlasreg"]["level"],
            suffix="spim.nii",
        ),
        xfm_ras=rules.affine_reg.output.xfm_ras,
    output:
        warp=bids(
            root=root,
            subject="{subject}",
            datatype="warps",
            sample="{sample}",
            acq="{acq}",
            from_="subject",
            to="{template}",
            suffix="warp.nii",
        ),
        warped=bids(
            root=root,
            subject="{subject}",
            datatype="warps",
            sample="{sample}",
            acq="{acq}",
            space="{template}",
            desc="deformwarped",
            suffix="spim.nii",
        ),
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="deform_reg",
            sample="{sample}",
            acq="{acq}",
            space="{template}",
            suffix="log.txt",
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
        zarr_zip=bids(
            root=root,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            stain=config["atlasreg"]["stain"],
            suffix="spim.ome.zarr.zip",
        ),
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="{desc}",
                    from_="{template}",
                    suffix="dseg.zarr",
                )
            )
        ),
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="resample_labels_to_zarr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            space="{template}",
            suffix="log.txt",
        ),
    script:
        "../scripts/resample_labels_to_zarr.py"


rule zarr_to_ome_zarr_labels:
    input:
        zarr=bids(
            root=work,
            subject="{subject}",
            datatype="micr",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            from_="{template}",
            suffix="dseg.zarr",
        ),
        metadata_json=rules.raw_to_metadata.output.metadata_json,
        label_tsv=bids_tpl(root=root, template="{template}", suffix="dseg.tsv"),
    params:
        max_downsampling_layers=config["ome_zarr"]["max_downsampling_layers"],
        rechunk_size=config["ome_zarr"]["rechunk_size"],
        scaling_method="nearest",
        downsampling=config["bigstitcher"]["fuse_dataset"]["downsampling"],
        label_name="dseg_{template}",
    output:
        zarr=temp(
            directory(
                bids(
                    root=work,
                    subject="{subject}",
                    datatype="micr",
                    sample="{sample}",
                    acq="{acq}",
                    desc="{desc}",
                    from_="{template}",
                    suffix="dseg.ome.zarr",
                )
            )
        ),
    threads: 32
    container:
        config["containers"]["spimprep"]
    group:
        "preproc"
    log:
        bids(
            root="logs",
            subject="{subject}",
            datatype="zarr_to_ome_zarr_labels",
            sample="{sample}",
            acq="{acq}",
            desc="{desc}",
            space="{template}",
            suffix="log.txt",
        ),
    script:
        "../scripts/zarr_to_ome_zarr_labels.py"
