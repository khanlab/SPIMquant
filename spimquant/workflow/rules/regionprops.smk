rule compute_filtered_regionprops:
    """Calculate region props from filtered objects of segmentation."""
    input:
        mask=bids_oz_in(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    output:
        regionprops_parquet=temp(
            bids(
                root=root,
                datatype="tabular",
                stain="{stain}",
                desc="{desc}",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=scaled_mem_mb(32000, 60, config["segmentation_level"]),
        runtime=scaled_runtime(30, 0.2, config["segmentation_level"]),
    params:
        region_filters=lambda wildcards: config["stain_defaults"][
            "regionprop_filters"
        ].get(wildcards.stain, config["regionprop_filters"]),
        output_properties=config["regionprop_outputs"],
    script:
        "../scripts/compute_filtered_regionprops.py"


rule transform_regionprops_to_template:
    """Transform regionprops coordinates from subject to template space."""
    input:
        regionprops_parquet=bids(
            root=root,
            datatype="tabular",
            stain="{stain}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        xfm_composite_inv=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to="subject",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        regionprops_transformed_parquet=temp(
            bids(
                root=root,
                datatype="tabular",
                stain="{stain}",
                desc="{desc}",
                space="{template}",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
    params:
        coord_column_names=config["coord_column_names"],
    script:
        "../scripts/transform_regionprops_to_template.py"


def get_regionprops_parquets_to_aggregate(wildcards):
    """Per-stain regionprops for the stains this method actually segmented.

    Cannot be resolved at parse time: a single-stain method such as the plaque
    ensemble aggregates over only its own stain, not all of stains_for_seg.
    """
    paths = expand(
        bids(
            root=root,
            datatype="tabular",
            stain="{stain}",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        stain=stains_for_desc(wildcards.desc),
        allow_missing=True,
    )
    return [path.format(**wildcards) for path in paths]


rule aggregate_regionprops_across_stains:
    """Aggregate transformed regionprops across stains."""
    input:
        regionprops_parquets=get_regionprops_parquets_to_aggregate,
    output:
        regionprops_aggregated_parquet=bids(
            root=root,
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    params:
        # must stay aligned with regionprops_parquets, which is per-method
        stains=lambda wildcards: stains_for_desc(wildcards.desc),
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
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    output:
        coloc_parquet=bids(
            root=root,
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="coloc.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 16
    resources:
        mem_mb=64000,
        runtime=180,
    params:
        coord_column_names=config["coord_column_names"],
        template_coord_column_names=config["template_coord_column_names"],
        search_radius_multiplier=1.0,
        overlap_threshold=0.0,
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
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        scalar=bids_oz_in(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level=config["segmentation_level"],
            desc=config["vessel_seg_method"],
            suffix="dist.{ext}",
            **inputs["spim"].wildcards,
        ),
    output:
        parquet=bids(
            root=root,
            datatype="tabular",
            desc="{desc}+vessels",
            space="{template}",
            vessels="{stain}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    threads: 8
    resources:
        mem_mb=4000,
        runtime=scaled_runtime(15, 0.05, config["segmentation_level"]),
    params:
        coord_column_names=config["coord_column_names"],
        col_name="vesseldist",
        zarrnii_kwargs=zarrnii_in_kwargs,
    script:
        "../scripts/sample_at_points.py"
