rule compute_filtered_regionprops:
    """Objects as connected components of the mask (--object_method components)."""
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
                obj="components",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        # dask gives each worker an equal share of this; with less, workers
        # pause at their limit and the job runs ~10x slower
        mem_mb=scaled_mem_mb(96000, 30, config["segmentation_level"]),
        runtime=scaled_runtime(30, 0.2, config["segmentation_level"]),
    params:
        region_filters=lambda wildcards: config["stain_defaults"][
            "regionprop_filters"
        ].get(wildcards.stain, config["regionprop_filters"]),
        output_properties=config["regionprop_outputs"],
    script:
        "../scripts/compute_filtered_regionprops.py"


rule compute_instance_regionprops:
    """Objects by seeded watershed of the mask (--object_method watershed).

    Splits the mask into instances with a seeded watershed: seeds are
    prominent local maxima of the smoothed corrected image, and every mask
    voxel goes to one seed, so touching puncta are counted separately
    instead of merging into one connected component.  Writes the same
    columns as compute_filtered_regionprops plus seed position, intensity
    and elongation, and a JSON file with the parameters used (including the
    seed prominence h and how it was set).
    """
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
        corrected=bids_oz_in(
            root=work,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="corrected{method}".format(method=config["correction_method"]),
            suffix="SPIM.{ext}",
            **inputs["spim"].wildcards,
        ),
    output:
        regionprops_parquet=temp(
            bids(
                root=root,
                datatype="tabular",
                stain="{stain}",
                desc="{desc}",
                obj="watershed",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
        stats_json=bids(
            root=root,
            datatype="tabular",
            stain="{stain}",
            desc="{desc}",
            obj="watershed",
            suffix="instances.json",
            **inputs["spim"].wildcards,
        ),
    threads: 64 if config["dask_scheduler"] == "distributed" else 32
    resources:
        # each task (one instance_chunk_size block plus overlap) peaks at
        # ~100 bytes per voxel of the overlapped block, mostly h_maxima. The
        # overlap is instance_depth_um on each side, so in voxels it grows as
        # the voxels shrink: at 1.5 x 0.975 x 0.975 um, 160^3 blocks become
        # 228 x 264 x 264 (~16M voxels), ~1.6 GB per task. With ~1.7 GB per
        # thread (128 threads, 16 workers of 13.4 GB), dense blocks pushed
        # workers past 95% and killed them repeatedly (KilledWorker) on a
        # 250 GB marmoset hemisphere. So the job uses 64 threads and ~3 GB per
        # thread (8 workers of 8 threads get ~24 GB each), whatever the image
        # size. What grows with the image is the scheduler's task graph
        # (~0.4 GB serialised for 31k blocks), reserved with the script's
        # memory.
        mem_mb=scaled_mem_mb(208000, 40, config["segmentation_level"]),
        runtime=scaled_runtime(60, 0.25, config["segmentation_level"]),
    params:
        dask_threads_per_worker=config["instance_threads_per_worker"],
        dask_reserve_mb=lambda wildcards: scaled_mem_mb(
            16000, 40, config["segmentation_level"]
        )(wildcards, 1),
        sigma_um=config["instance_sigma_um"],
        min_dist_um=config["instance_min_dist_um"],
        h=config["instance_h"],
        noise_k=config["instance_h_noise_k"],
        noise_percentile=config["instance_h_noise_percentile"],
        h_fallback=config["instance_h_fallback"],
        bg_sigma_um=config["instance_bg_sigma_um"],
        stat_blocks=config["instance_stat_blocks"],
        stat_samples_per_block=config["instance_stat_samples_per_block"],
        stat_level_offset=config["instance_stat_level_offset"],
        min_volume_um3=config["instance_min_volume_um3"],
        depth_um=config["instance_depth_um"],
        chunk_size=config["instance_chunk_size"],
        task_batch=config["instance_task_batch"],
        output_properties=config["regionprop_outputs"],
    script:
        "../scripts/compute_instance_regionprops.py"


rule transform_regionprops_to_template:
    """Add template-space coordinates to a stain's objects.

    Also labels each row with its stain and methods (``stain``,
    ``seg_method``, ``object_method``), so per-stain and aggregated tables
    have the same columns.
    """
    input:
        regionprops_parquet=bids(
            root=root,
            datatype="tabular",
            stain="{stain}",
            desc="{desc}",
            obj="{obj}",
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
                obj="{obj}",
                space="{template}",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj_quant"
    threads: 1
    resources:
        mem_mb=scaled_mem_mb(2000, 18000, config["registration_level"]),
        runtime=10,
    params:
        coord_column_names=config["coord_column_names"],
        labels=lambda wildcards: {
            "stain": wildcards.stain,
            "seg_method": wildcards.desc,
            "object_method": wildcards.obj,
        },
    script:
        "../scripts/transform_regionprops_to_template.py"


def primary_regionprops(stain, template_space=True):
    """Path pattern of a stain's objects from its primary methods."""
    kwargs = {"space": "{template}"} if template_space else {}
    return bids(
        root=root,
        datatype="tabular",
        stain=stain,
        desc=seg_plan[stain].seg_method,
        obj=seg_plan[stain].object_method,
        suffix="regionprops.parquet",
        **kwargs,
        **inputs["spim"].wildcards,
    )


rule aggregate_regionprops_across_stains:
    """Objects of every stain with objects, from each stain's primary methods."""
    input:
        regionprops_parquets=[primary_regionprops(s) for s in stains_with_objects],
    output:
        regionprops_aggregated_parquet=temp(
            bids(
                root=root,
                datatype="tabular",
                space="{template}",
                desc="merged",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj_quant"
    threads: 1
    resources:
        mem_mb=1500,
        runtime=10,
    script:
        "../scripts/aggregate_regionprops_across_stains.py"


# stains whose mask's signed distance transform is sampled at every object
# (<stain>_sdt): the vessel stain, and any of the sdt_stains config key that
# are segmented (from their primary mask)
sdt_stains = ([stain_for_vessels] if do_vessels else []) + [
    s
    for s in (config.get("sdt_stains") or [])
    if s in seg_plan and s != stain_for_vessels
]


def sdt_input(wildcards):
    """Signed distance transform of ``wildcards.sdt_stain``'s mask."""
    stain = wildcards.sdt_stain
    if do_vessels and stain == stain_for_vessels:
        datatype, desc = "vessels", config["vessel_seg_method"]
    else:
        datatype, desc = "seg", seg_plan[stain].seg_method
    return bids_oz_in(
        root=root,
        datatype=datatype,
        stain=stain,
        level=config["segmentation_level"],
        desc=desc,
        suffix="dist.{ext}",
        **inputs["spim"].wildcards,
    ).format(**wildcards)


rule sample_sdt_at_objects:
    """Signed distance to a stain's mask at every object centroid (mm).

    Negative inside the mask, positive outside. Writes only the new column
    (<stain>_sdt), in the row order of the aggregated objects.
    """
    input:
        parquet=rules.aggregate_regionprops_across_stains.output.regionprops_aggregated_parquet,
        scalar=sdt_input,
    output:
        parquet=temp(
            bids(
                root=root,
                datatype="tabular",
                space="{template}",
                desc="{sdt_stain}sdt",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            )
        ),
    wildcard_constraints:
        sdt_stain="|".join(sdt_stains) or "^$",
    group:
        "subj_quant"
    threads: 8
    resources:
        mem_mb=4000,
        runtime=scaled_runtime(15, 0.05, config["segmentation_level"]),
    params:
        coord_column_names=config["coord_column_names"],
        col_name=lambda wildcards: f"{wildcards.sdt_stain}_sdt",
        column_only=True,
        # the SDT is written by SPIMquant; the raw-input overrides (orientation,
        # --sample_staining channel labels) do not apply to it
        zarrnii_kwargs={},
    script:
        "../scripts/sample_at_points.py"


rule annotate_object_distances:
    """The objects of all stains, with distances between them (all in mm).

    For every stain S with objects, each object gets S_dist, the distance from
    its centroid to the nearest centroid of S (the nearest other object for
    its own stain), and S_edgedist, that minus the nearest object's
    equivalent-sphere radius (from its volume), clipped at 0: an approximate
    distance to its surface. The sampled signed distances (<stain>_sdt) are
    added too.
    """
    input:
        regionprops_parquet=rules.aggregate_regionprops_across_stains.output.regionprops_aggregated_parquet,
        sdt_parquets=expand(
            bids(
                root=root,
                datatype="tabular",
                space="{template}",
                desc="{sdt_stain}sdt",
                suffix="regionprops.parquet",
                **inputs["spim"].wildcards,
            ),
            sdt_stain=sdt_stains,
            allow_missing=True,
        ),
    output:
        regionprops_parquet=bids(
            root=root,
            datatype="tabular",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_quant"
    threads: 4
    resources:
        mem_mb=scaled_mem_mb(4000, 30, config["segmentation_level"]),
        runtime=15,
    params:
        coord_column_names=config["coord_column_names"],
    script:
        "../scripts/object_distances.py"


rule colocalize_regionprops:
    """Colocalization across every pair of stains with objects.

    Optional parameters (with defaults in script):
        - search_radius_multiplier: 1.0 (controls search distance)
        - overlap_threshold: 0.0 (minimum overlap to record)
    """
    input:
        regionprops_parquet=rules.aggregate_regionprops_across_stains.output.regionprops_aggregated_parquet,
    output:
        coloc_parquet=bids(
            root=root,
            datatype="tabular",
            space="{template}",
            suffix="coloc.parquet",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_coloc"
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
