
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



rule coiled_applyn4_otsu:
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
        cluster_name=bids(
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                suffix="n4_otsu",
                **inputs["spim"].wildcards
            ),
        otsu_n_classes=3,
        otsu_threshold_index=-1, #-1 selects the highest intensity threshold from all the classes
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
        mask_uri=bids(
                root=root_coiled,
                datatype="micr",
                stain="{stain}",
                dslevel="{dslevel}",
                level="{level}",
                desc="n4otsu",
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
                desc="n4otsu",
                suffix="mask.DONE",
                **inputs["spim"].wildcards))
    threads: 1
    container: None
    script: '../scripts/coiled_n4_otsu.py'



