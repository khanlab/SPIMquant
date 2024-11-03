
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
        cluster_name='n4-{subject}',
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
    container: None
    script: '../scripts/coiled_n4.py'


rule coiled_otsu:
    input:
        rules.coiled_n4.output
    params:
        cluster_name='otsu-{subject}',
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
    script:
        "../scripts/coiled_fieldfrac.py"
