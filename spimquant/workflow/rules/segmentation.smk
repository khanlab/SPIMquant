
rule gaussian_biasfield:
    """simple bias field correction with gaussian"""
    input:
        spim=inputs["spim"].path,
    params:
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        corrected=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="corrected",
                corrmethod="gaussian",
                suffix="SPIM.ome.zarr",
                **inputs["spim"].wildcards,
            )
        ),
        biasfield=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="gaussian",
                suffix="biasfield.ome.zarr",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 128
    script:
        "../scripts/gaussian_biasfield.py"


rule n4_biasfield:
    """N4 bias field correction with antspyx"""
    input:
        spim=inputs["spim"].path,
    params:
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        corrected=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="corrected",
                corrmethod="n4",
                suffix="SPIM.ome.zarr",
                **inputs["spim"].wildcards,
            )
        ),
        biasfield=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="n4",
                suffix="biasfield.ome.zarr",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 128
    script:
        "../scripts/n4_biasfield.py"


rule multiotsu:
    input:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
            level="{level}",
            desc="corrected",
            corrmethod=config["correction_method"],
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
        mask=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="otsu+k{k,[0-9]+}i{i,[0-9]+}",
                suffix="mask.ome.zarr",
                **inputs["spim"].wildcards,
            )
        ),
        thresholds_png=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
            level="{level}",
            desc="otsu+k{k,[0-9]+}i{i,[0-9]+}",
            suffix="thresholds.png",
            **inputs["spim"].wildcards,
        ),
    threads: 128
    script:
        "../scripts/multiotsu.py"


rule threshold:
    input:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
            level="{level}",
            desc="corrected",
            corrmethod=config["correction_method"],
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
                dslevel="{dslevel}",
                level="{level}",
                desc="threshold",
                suffix="mask.DONE",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 128
    script:
        "../scripts/threshold.py"


rule clean_segmentation:
    input:
        mask=bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="{desc}",
                suffix="mask.ome.zarr",
                **inputs["spim"].wildcards,
        ),
    params:
        max_extent=0.15,
        zarrnii_kwargs={}
    output:
        exclude_mask=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="{desc}+cleaned",
                suffix="excludemask.ome.zarr",
                **inputs["spim"].wildcards,
            )
        ),
        cleaned_mask=directory(bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="{desc}+cleaned",
                suffix="mask.ome.zarr",
                **inputs["spim"].wildcards,
                )
        ),
    threads: 128
    script:
        "../scripts/clean_segmentation.py"

   
rule compute_centroids:
    """Calculate object/cell centroids from the segmentation"""
    input:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel=config["registration_level"],
            level=config["segmentation_level"],
            desc="{seg_method}",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        zarrnii_kwargs={}
    output:
        centroids_npy=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
            desc="{seg_method}",
            suffix="centroids.npy",
            **inputs["spim"].wildcards,
        ),
    threads: 128
    script:
        "../scripts/compute_centroids.py"



rule fieldfrac:
    """ Calculates fieldfrac from a binary mask via downsampling, assuming mask intensity is 100"""
    input:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel=config["registration_level"],
            level=config["segmentation_level"],
            desc="{seg_method}",
            suffix="mask.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    params:
        ds_level=lambda wildcards: int(wildcards.dslevel)
        - int(config["segmentation_level"]),
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        fieldfrac_nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
            desc="{seg_method}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards,
        ),
    threads: 32
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
    threads: 32
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        " -ri NN "
        "  -rm {input.mask} {output.mask} "
        "  -r {input.xfm_ras},-1 {input.invwarp}"
        #note: LABEL interpolation not possible with >1000 labels


rule map_img_to_roi_tsv:
    input:
        img=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
            desc="{desc}",
            suffix="{suffix}.nii",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{dslevel}",
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
            dslevel="{dslevel}",
            desc="{desc}",
            suffix="{suffix,fieldfrac}stats.tsv",
            **inputs["spim"].wildcards,
        ),
    script:
        "../scripts/map_img_to_roi_tsv.py"


rule map_centroids_to_atlas_rois:
    input:
        centroids_npy=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
            desc="{seg_method}",
            suffix="centroids.npy",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{dslevel}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    output:
        centroids_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            dslevel="{dslevel}",
            desc="{seg_method}",
            suffix="centroidstats.tsv",
            **inputs["spim"].wildcards,
        ),
        counts_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            dslevel="{dslevel}",
            desc="{seg_method}",
            suffix="countstats.tsv",
            **inputs["spim"].wildcards,
        ),

    script:
        "../scripts/map_atlas_to_centroids.py"


rule merge_into_segstats_tsv:
    input:
        counts_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            dslevel="{dslevel}",
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
            dslevel="{dslevel}",
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
            dslevel="{dslevel}",
            desc="{desc}",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards,
        ),
    script: '../scripts/merge_into_segstats_tsv.py'

rule map_segstats_tsv_dseg_to_template_nii:
    """ uses generic script that paints regions with column data (e.g. use this to make density heat-maps)"""
    input:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            dslevel=config["registration_level"],
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
            dslevel="{dslevel}",
            desc="{desc}",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{dslevel}",
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
            dslevel="{dslevel}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="{suffix}.nii",
            **inputs["spim"].wildcards,
        ),
    script:
        "../scripts/map_tsv_dseg_to_nii.py"


rule deform_fieldfrac_nii_to_template_nii:
    input:
        flo=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            dslevel="{dslevel}",
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
            dslevel="{dslevel}",
            desc="{desc}",
            space="{template}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.flo} -o {output.nii} "
        " -r {input.ref} -t {input.warp} {input.xfm_itk}"
