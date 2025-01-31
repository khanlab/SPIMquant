
rule antspyx_n4:
    input:
        **get_storage_creds(inputs["spim"].path),
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



rule dask_n4:
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
    threads: 1 if config['use_coiled'] else 32
    resources: 
        coiled=1 
    container: None
    script: '../scripts/dask_n4.py'


rule downsampled_apply_n4_mask:
    input:
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
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config['registration_level'],
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
        "c3d {input.n4_bf_ds} {input.spim_ds} -divide -as N4 -replace inf 10000  {input.mask} -reslice-identity -push N4 -multiply -o {output.masked}"

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
    params:
        otsu_n_classes=3,
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
    script:
        '../scripts/calc_otsu_thresholds.py'


rule dask_otsu:
    input:
        n4=rules.dask_n4.output,
        otsu_thresholds=bids(
                root=work,
                datatype="micr",
                stain="{stain}",
                level="{dslevel}",
                desc="n4corrmasked",
                suffix="thresholds.npy",
                **inputs["spim"].wildcards
            ),
    params:
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
    threads: 1 if config['use_coiled'] else 32
    resources: 
        coiled=1 
    container: None
    script: '../scripts/dask_otsu.py'



rule dask_fieldfrac:
    input:
        bids(root=root,
                datatype="micr",
                stain="{stain}",
                dslevel=config['registration_level'],
                level=config['segmentation_level'],
                desc="otsu",
                suffix="mask.DONE",
                **inputs["spim"].wildcards)
    params:
        mask_uri=bids(
                root=root_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel=config['registration_level'],
                level=config['segmentation_level'],
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
        "../scripts/dask_fieldfrac.py"


rule deform_negative_mask_to_subject_nii:
    input:
        ref=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            suffix="SPIM.nii",
            **inputs["spim"].wildcards
        ),
        mask=config['template_negative_mask'],
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
            **inputs["spim"].wildcards
        ),
    threads: 32
    container:
        config["containers"]["itksnap"]
    shell:
        " greedy -threads {threads} -d 3 -rf {input.ref} "
        " -ri NN "
        "  -rm {input.mask} {output.mask} "
        "  -r {input.xfm_ras},-1 {input.invwarp}"
        #note: LABEL interpolation not possible with >1000 labels


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
        mask1=bids(
            root=root,
            datatype="micr",
            desc="negative",
            level="{dslevel}",
            from_=config["template"],
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards
        ),
        mask2=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            desc="brain",
            level="{dslevel}",
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
        "c3d {input.fieldfrac} {input.mask1} -multiply {input.mask2} -multiply -o {output.fieldfrac_mod}"

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
                desc="otsupenalty" if config['use_negative_mask'] else 'otsu',
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
            desc="otsupenalty" if config['use_negative_mask'] else 'otsu',
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
            dslevel=config["registration_level"],
            desc="otsupenalty" if config['use_negative_mask'] else 'otsu',
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
            desc="otsupenalty" if config['use_negative_mask'] else 'otsu',
            suffix="segstats.tsv",
            **inputs["spim"].wildcards
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{dslevel}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards
        ),
    params:
        label_column="index",
        feature_column="avg_fieldfrac",
    output:
        nii=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            dslevel="{dslevel}",
            from_="{template}",
            stain="{stain}",
            suffix="fieldfrac.nii",
            **inputs["spim"].wildcards
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
                **inputs["spim"].wildcards),
        ref=bids_tpl(root=root, template="{template}", desc="LR", suffix="dseg.nii.gz"),
        xfm_itk=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards
        ),
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii",
            **inputs["spim"].wildcards
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
                **inputs["spim"].wildcards),
    threads: 32
    container:
        config["containers"]["ants"]
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.flo} -o {output.nii} "
        " -r {input.ref} -t {input.warp} {input.xfm_itk}"



