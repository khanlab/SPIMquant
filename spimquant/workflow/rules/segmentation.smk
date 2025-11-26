
rule gaussian_biasfield:
    """simple bias field correction with gaussian"""
    input:
        spim=inputs["spim"].path,
    params:
        proc_level=5,
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        corrected=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="correctedgaussian",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
        biasfield=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="gaussian",
                    suffix="biasfield.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=15,
    script:
        "../scripts/gaussian_biasfield.py"


rule n4_biasfield:
    """N4 bias field correction with antspyx"""
    input:
        spim=inputs["spim"].path,
    params:
        proc_level=5,
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        corrected=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="correctedn4",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
        biasfield=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="n4",
                    suffix="biasfield.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=500000,
        disk_mb=2097152,
        runtime=60,
    script:
        "../scripts/n4_biasfield.py"


rule multiotsu:
    input:
        corrected=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="corrected{method}".format(method=config["correction_method"]),
            suffix="SPIM.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        hist_bins=int(config["seg_hist_bins"]),
        hist_range=[int(x) for x in config["seg_hist_range"]],
        otsu_k=lambda wildcards: int(wildcards.k),
        otsu_threshold_index=lambda wildcards: int(wildcards.i),
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        mask=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="otsu+k{k,[0-9]+}i{i,[0-9]+}",
                    suffix="mask.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
        thresholds_png=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="otsu+k{k,[0-9]+}i{i,[0-9]+}",
            suffix="thresholds.png",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=500000,
        disk_mb=2097152,
        runtime=15,
    script:
        "../scripts/multiotsu.py"


rule convert_zarr_to_ozx:
    """generic rule to convert ome zarr to zip (.ozx)"""
    input:
        zarr=str(Path(work) / "{prefix}.ome.zarr"),
    output:
        ozx=str(Path(root) / "{prefix}.ozx"),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=60,
    group:
        "subj"
    script:
        "../scripts/convert_zarr_to_ozx.py"


rule threshold:
    input:
        corrected=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="corrected{method}".format(method=config["correction_method"]),
            suffix="SPIM.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        threshold=int(config["seg_threshold"]),
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        mask=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="threshold",
                suffix="mask.DONE",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=256000,
        runtime=15,
    script:
        "../scripts/threshold.py"


rule clean_segmentation:
    input:
        mask=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        max_extent=0.15,
        proc_level=2,  #level at which to calculate conncomp
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        exclude_mask=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="{desc}+cleaned",
                    suffix="excludemask.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
        cleaned_mask=temp(
            directory(
                bids(
                    root=work,
                    datatype="micr",
                    stain="{stain}",
                    level="{level}",
                    desc="{desc}+cleaned",
                    suffix="mask.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            )
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=30,
    script:
        "../scripts/clean_segmentation.py"


rule compute_filtered_regionprops:
    """Calculate region props from filtered objects of segmentation."""
    input:
        mask=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        region_filters=config["regionprop_filters"],
        output_properties=config["regionprop_outputs"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        regionprops_parquet=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 128
    resources:
        mem_mb=256000,
        runtime=30,
    script:
        "../scripts/compute_filtered_regionprops.py"


rule counts_per_voxel:
    """Calculate counts per voxel based on points"""
    input:
        ref_spim=inputs["spim"].path,
        regionprops_parquet=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    params:
        coord_column_names=config["coord_column_names"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        counts_nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="counts.nii",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 16
    resources:
        mem_mb=15000,
        runtime=10,
    script:
        "../scripts/counts_per_voxel.py"


rule fieldfrac:
    """ Calculates fieldfrac from a binary mask via downsampling, assuming mask intensity is 100.
        The level in the input corresponds to the level of the input mask, and the level in the 
        output image is the level of the downsampled fieldfrac map. Internally we calculate 
        what downsampling factor to use to achieve the desired level
        """
    input:
        mask=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        hires_level=config["segmentation_level"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        fieldfrac_nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/fieldfrac.py"


rule deform_negative_mask_to_subject_nii:
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards,
        ),
        mask=config["template_negative_mask"],
        xfm_ras=rules.init_affine_reg.output.xfm_ras,
        invwarp=rules.deform_reg.output.invwarp,
    output:
        mask=bids(
            root=root,
            datatype="micr",
            desc="negative",
            level="{level}",
            from_="{template}",
            suffix="mask.nii.gz",
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
        " -ri NN "
        "  -rm {input.mask} {output.mask} "
        "  -r {input.xfm_ras},-1 {input.invwarp}"


rule map_img_to_roi_tsv:
    input:
        img=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="{suffix}.nii",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{level}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    output:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="{suffix,fieldfrac}stats.tsv",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/map_img_to_roi_tsv.py"


rule map_regionprops_to_atlas_rois:
    input:
        regionprops_parquet=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{level}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        coord_column_names=config["coord_column_names"],
    output:
        regionprops_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="regionpropstats.tsv",
            **inputs["spim"].wildcards,
        ),
        counts_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="countstats.tsv",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/map_atlas_to_regionprops.py"


rule merge_into_segstats_tsv:
    input:
        regionprops_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="regionpropstats.tsv",
            **inputs["spim"].wildcards,
        ),
        counts_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="countstats.tsv",
            **inputs["spim"].wildcards,
        ),
        fieldfrac_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="fieldfracstats.tsv",
            **inputs["spim"].wildcards,
        ),
    output:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/merge_into_segstats_tsv.py"


rule map_segstats_tsv_dseg_to_template_nii:
    """ uses generic script that paints regions with column data (e.g. use this to make density heat-maps)"""
    input:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level=config["registration_level"],
            desc="{desc}",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards,
        ),
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        label_column="index",
        feature_column="{suffix}",
    output:
        nii=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            space="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="{suffix}.nii",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/map_tsv_dseg_to_nii.py"


rule map_segstats_tsv_dseg_to_subject_nii:
    """ uses generic script that paints regions with column data (e.g. use this to make density heat-maps)"""
    input:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{level}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        label_column="index",
        feature_column="{suffix}",
    output:
        nii=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            level="{level}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="{suffix}.nii",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/map_tsv_dseg_to_nii.py"


rule deform_fieldfrac_nii_to_template_nii:
    input:
        flo=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards,
        ),
        ref=rules.import_template_anat.output.anat,
        xfm_itk=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii",
            **inputs["spim"].wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            space="{template}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.flo} -o {output.nii} "
        " -r {input.ref} -t {input.warp} {input.xfm_itk}"
