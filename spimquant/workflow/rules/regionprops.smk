rule compute_filtered_regionprops:
    """Calculate region props from filtered objects of segmentation."""
    input:
        mask=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    params:
        region_filters=lambda wildcards: config["stain_defaults"][
            "regionprop_filters"
        ].get(wildcards.stain, config["regionprop_filters"]),
        output_properties=config["regionprop_outputs"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        regionprops_parquet=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                desc="{desc}",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=256000,
        runtime=180,
    script:
        "../scripts/compute_filtered_regionprops.py"


rule transform_regionprops_to_template:
    """Transform regionprops coordinates from subject to template space."""
    input:
        regionprops_parquet=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        xfm_ras=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="ras",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
        invwarp=bids(
            root=root,
            datatype="warps",
            from_="{template}",
            to="subject",
            suffix="warp.nii.gz",
            **inputs["spim"].wildcards,
        ),
    params:
        coord_column_names=config["coord_column_names"],
    output:
        regionprops_transformed_parquet=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                desc="{desc}",
                space="{template}",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/transform_regionprops_to_template.py"


rule aggregate_regionprops_across_stains:
    """Aggregate transformed regionprops across stains."""
    input:
        regionprops_parquets=expand(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                desc="{desc}",
                space="{template}",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            ),
            stain=stains_for_seg,
            allow_missing=True,
        ),
    params:
        stains=stains_for_seg,
    output:
        regionprops_aggregated_parquet=bids(
            root=root,
            datatype="micr",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    script:
        "../scripts/aggregate_regionprops_across_stains.py"


rule colocalize_regionprops:
    """Perform colocalization analysis across channel pairs.
    
    Optional parameters (with defaults in script):
        - search_radius_multiplier: 1.0 (controls search distance)
        - overlap_threshold: 0.0 (minimum overlap to record)
    """
    input:
        regionprops_parquet=bids(
            root=root,
            datatype="micr",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    params:
        coord_column_names=config["coord_column_names"],
        template_coord_column_names=config["template_coord_column_names"],
        search_radius_multiplier=1.0,
        overlap_threshold=0.0,
    output:
        coloc_parquet=bids(
            root=root,
            datatype="micr",
            desc="{desc}",
            space="{template}",
            suffix="coloc.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=30,
    script:
        "../scripts/compute_colocalization.py"


rule sample_at_vessel_sdt:
    """Sample sdt at points
    - < 0 : instance is inside the mask
    - > 0 : instance is outside the mask

    """
    input:
        parquet=bids(
            root=root,
            datatype="micr",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        scalar=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=config["segmentation_level"],
            desc=config["vessel_seg_method"],
            suffix="dist.ozx",
            **inputs["spim"].wildcards,
        ),
    params:
        coord_column_names=config["coord_column_names"],
        col_name="sdt_{stain}",
        zarrnii_kwargs={"orientation": config["orientation"], "level": 0},
    output:
        parquet=bids(
            root=root,
            datatype="micr",
            desc="{desc}",
            space="{template}",
            vessels="{stain}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=32000,
        runtime=30,
    script:
        "../scripts/sample_at_points.py"


