
rule antspyx_n4:
    params:
        spim_uri=inputs["spim"].path,
    output:
        spim_ds=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                suffix="SPIM.nii",
                **inputs["spim"].wildcards
            ),
        n4_bf_ds=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                suffix="n4biasfield.nii",
                **inputs["spim"].wildcards
            ),
    threads: 8
    resources:
        mem_mb=16000
    container: None
    script: '../scripts/antspyx_n4.py'



rule coiled_n4:
    input:
        n4_bf_ds=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{dslevel}",
                suffix="n4biasfield.nii",
                **inputs["spim"].wildcards
            ),
    params:
        spim_uri=inputs["spim"].path,
        bf_ds_uri=bids(
                root=work_coiled,
                datatype="micr",
                stain="{stain}",
                level="{dslevel}",
                suffix="n4biasfield.ome.zarr",
                **inputs["spim"].wildcards
            ),
        bf_us_uri=bids(
                root=work_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                suffix="n4biasfield.ome.zarr",
                **inputs["spim"].wildcards
            ),
        spim_n4_uri=bids(
                root=root_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="n4corr",
                suffix="SPIM.ome.zarr",
                **inputs["spim"].wildcards
            ),
    output:
            touch(bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="n4corr",
                suffix="SPIM.DONE",
                **inputs["spim"].wildcards))
    threads: 1
    resources: 
        coiled=1 
    container: None
    script: '../scripts/coiled_n4.py'


rule get_downsampled_n4:
    input:
        coiled_n4=bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel=config["segment"]["n4_ds_level"],
                level=0,
                desc="n4corr",
                suffix="SPIM.DONE",
                **inputs["spim"].wildcards),

    params:
        spim_n4_uri=bids(
                root=root_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel=config["segment"]["n4_ds_level"],
                level=0,
                desc="n4corr",
                suffix="SPIM.ome.zarr",
                **inputs["spim"].wildcards
            ),
    output:
        nii=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="n4corr",
                suffix="SPIM.nii",
                **inputs["spim"].wildcards
            ),

rule mask_downsampled_n4:
    input:
        corrected=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="n4corr",
                suffix="SPIM.nii",
                **inputs["spim"].wildcards
            ),
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="brain",
            suffix="mask.nii",
            **inputs["spim"].wildcards
        ),
    output:
        masked=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="n4corrmasked",
                suffix="SPIM.nii",
                **inputs["spim"].wildcards
            ),

    container:
        config["containers"]["itksnap"]
    shell:
        "c3d {input.corrected} {input.mask} -multiply -o {output.masked}"

#calc thresholds using downsampled, masked image
rule calc_otsu_thresholds:
    input:
        masked=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="n4corrmasked",
                suffix="SPIM.nii",
                **inputs["spim"].wildcards
            ),
    output:
        otsu_thresholds=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{level}",
                desc="n4corrmasked",
                suffix="thresholds.npy",
                **inputs["spim"].wildcards
            ),



#TODO: try with fixed threshold, as some images pre-processing seems to have caused issues (ie fieldfrac correction failed)
rule coiled_otsu:
    input:
        rules.coiled_n4.output
    params:
        otsu_n_classes=3,
        otsu_threshold_index=-1, #-1 selects the highest intensity threshold from all the classes
        spim_n4_uri=bids(
                root=root_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="n4corr",
                suffix="SPIM.ome.zarr",
                **inputs["spim"].wildcards
            ),
        mask_uri=bids(
                root=root_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="otsu",
                suffix="mask.ome.zarr",
                **inputs["spim"].wildcards
            )
    output:
            touch(bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="otsu",
                suffix="mask.DONE",
                **inputs["spim"].wildcards))
    threads: 1
    resources: 
        coiled=1 
    container: None
    script: '../scripts/coiled_otsu.py'



rule coiled_fieldfrac:
    input:
        bids(root=root,
                datatype="micr",
                stain="{stain}",
                dslevel=config['segment']['n4_ds_level'],
                level=config['segment']['otsu_level'],
                desc="otsu",
                suffix="mask.DONE",
                **inputs["spim"].wildcards)
    params:
        mask_uri=bids(
                root=root_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel=config['segment']['n4_ds_level'],
                level=config['segment']['otsu_level'],
                desc="otsu",
                suffix="mask.ome.zarr",
                **inputs["spim"].wildcards
            )
    output:
        fieldfrac_nii=bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                desc="otsu",
                suffix="fieldfrac.nii",
                **inputs["spim"].wildcards)
    container: None
    threads: 1
    resources: 
        coiled=1 
    script:
        "../scripts/coiled_fieldfrac.py"


rule apply_boundary_penalty:
    input:
        fieldfrac=bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                desc="otsu",
                suffix="fieldfrac.nii",
                **inputs["spim"].wildcards),
        penalty=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{dslevel}",
            desc="brain",
            suffix="penalty.nii",
            **inputs["spim"].wildcards
        ),
    output:
        fieldfrac_mod=bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                desc="otsupenalty",
                suffix="fieldfrac.nii",
                **inputs["spim"].wildcards)
    container:
        config["containers"]["itksnap"]
    shell:
        "c3d {input.fieldfrac} -as FIELDFRAC {input.penalty} -reslice-identity -push FIELDFRAC -multiply -o {output.fieldfrac_mod}"

#now we have fieldfrac modulated by brainmask boundary penalty
#just need to calc avg fieldfrac in each ROI 
# if we then want total volume of plaques in each ROI, it is avg_fieldfrac * volume of voxel * number of voxels


rule map_fieldfrac_to_atlas_rois:
    input:
        img=bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                desc="otsupenalty",
                suffix="fieldfrac.nii",
                **inputs["spim"].wildcards),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level='{dslevel}',
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards
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
            desc="otsupenalty",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards
        ),
    script:
        "../scripts/map_img_to_roi_tsv.py"

rule map_segstats_tsv_dseg_to_template_nii:
    """ uses generic script that paints regions with column data (e.g. use this to make density heat-maps)"""
    input:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            dslevel=config["segment"]["fieldfrac_ds_level"],
            desc="otsupenalty",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards
        ),
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
    params:
        label_column="index",
        feature_column="avg_fieldfrac",
    output:
        nii=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            space="{template}",
            stain="{stain}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards
        ),
    script:
        "../scripts/map_tsv_dseg_to_nii.py"

