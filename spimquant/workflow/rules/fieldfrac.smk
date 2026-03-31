
rule fieldfrac:
    """Calculate field fraction from binary mask.
    
    Computes the fraction of brain tissue occupied by the segmented pathology at each
    voxel by downsampling the high-resolution mask. The output resolution (level) can
    differ from the input mask resolution, with the downsampling factor calculated
    automatically. Field fraction values range from 0-100.
    """
    input:
        mask=bids(
            root=root,
            datatype="{datatype}",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    params:
        hires_level=config["segmentation_level"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        fieldfrac_nii=bids(
            root=root,
            datatype="{datatype,seg|vessels}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="fieldfrac.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=30,
    script:
        "../scripts/fieldfrac.py"


rule map_fieldfrac_img_to_seg_tsv:
    input:
        img=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="{suffix}.nii.gz",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level="{level}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    output:
        tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                stain="{stain}",
                level="{level}",
                desc="{desc}",
                suffix="{suffix,fieldfrac}stats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    script:
        "../scripts/map_img_to_roi_tsv.py"
