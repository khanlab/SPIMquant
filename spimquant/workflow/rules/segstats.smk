"""Per-region statistics of the segmentations.

Per-stain tables carry the stain's methods in their names (desc-<seg method>,
obj-<object method>) and exist for every method of the stain, primary and
comparison. The merged table (mergedsegstats) has no method in its name: it
holds every segmented stain from its primary methods, the colocalization and
fibre ROI metrics, with ``<stain>+<metric>`` columns.
"""


def regionprops_for_roi_stats(wildcards):
    """The objects whose per-region statistics are wanted.

    Primary methods: the aggregated objects, which carry the distance columns
    (<stain>_dist, <stain>_edgedist, <stain>_sdt), so their regional means come
    through. Comparison methods: the stain's own objects.
    """
    if is_primary(wildcards.stain, wildcards.desc, wildcards.obj):
        return rules.annotate_object_distances.output.regionprops_parquet.format(
            **wildcards
        )
    return rules.transform_regionprops_to_template.output.regionprops_transformed_parquet.format(
        **wildcards
    )


rule map_regionprops_to_atlas_rois:
    input:
        regionprops_parquet=regionprops_for_roi_stats,
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
        regionprops_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            obj="{obj}",
            suffix="regionpropstats.tsv",
            **inputs["spim"].wildcards,
        ),
        counts_tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                stain="{stain}",
                level="{level}",
                desc="{desc}",
                obj="{obj}",
                suffix="countstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj_quant"
    threads: 4
    resources:
        mem_mb=scaled_mem_mb(2000, 10000),
        runtime=10,
    params:
        coord_column_names=config["coord_column_names"],
    script:
        "../scripts/map_atlas_to_regionprops.py"


rule map_coloc_to_atlas_rois:
    input:
        coloc_parquet=rules.colocalize_regionprops.output.coloc_parquet,
        dseg=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    output:
        coloc_tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                suffix="colocstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
        counts_tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                suffix="coloccountstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj_coloc"
    threads: 4
    resources:
        mem_mb=32000,
        runtime=15,
    params:
        coord_column_names=["template_coloc_x", "template_coloc_y", "template_coloc_z"],
    script:
        "../scripts/map_atlas_to_coloc.py"


def fieldfracstats_path(stain, desc, **kwargs):
    """Per-region field fraction table of ``stain``'s ``desc`` mask."""
    return bids(
        root=root,
        datatype="tabular",
        seg="{seg}",
        from_="{template}",
        stain=stain,
        level=config["registration_level"],
        desc=desc,
        suffix="fieldfracstats.tsv",
        **inputs["spim"].wildcards,
    )


rule merge_into_segstats_tsv:
    """Per-region statistics of one stain and one (seg method, object method)."""
    input:
        regionprops_tsv=rules.map_regionprops_to_atlas_rois.output.regionprops_tsv,
        counts_tsv=rules.map_regionprops_to_atlas_rois.output.counts_tsv,
        fieldfrac_tsv=bids(
            root=root,
            datatype="tabular",
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
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            obj="{obj}",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_quant"
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    script:
        "../scripts/merge_into_segstats_tsv.py"


def get_coloc_reference_fieldfrac(wildcards):
    """Fieldfrac table used only to recover ROI volume, so count becomes density.

    Every stain's field fraction table carries the same ROI volume column, so
    the first stain with objects is an arbitrary but stable choice.
    """
    stain = stains_with_objects[0]
    return fieldfracstats_path(stain, seg_plan[stain].seg_method).format(**wildcards)


rule merge_into_colocsegstats_tsv:
    """also includes fieldfracstats.tsv to obtain the volume to turn count into density"""
    input:
        coloc_tsv=rules.map_coloc_to_atlas_rois.output.coloc_tsv,
        counts_tsv=rules.map_coloc_to_atlas_rois.output.counts_tsv,
        fieldfrac_tsv=get_coloc_reference_fieldfrac,
    output:
        tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                suffix="colocsegstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj_coloc"
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    params:
        columns_to_drop=["fieldfrac"],
    script:
        "../scripts/merge_into_segstats_tsv.py"


def merged_table_sources(wildcards):
    """(prefix, table) pairs merged into mergedsegstats, in column order.

    Every segmented stain from its primary methods (its per-region statistics,
    or for --object_method none just its field fraction), then the fibre ROI
    statistics of each fibre stain (--fibres).
    """
    out = []
    for stain in stains_for_seg:
        p = seg_plan[stain]
        if p.has_objects:
            path = rules.merge_into_segstats_tsv.output.tsv.format(
                stain=stain,
                level=config["registration_level"],
                desc=p.seg_method,
                obj=p.object_method,
                **wildcards,
            )
        else:
            path = fieldfracstats_path(stain, p.seg_method).format(**wildcards)
        out.append((stain, path))
    for stain in stains_for_fibres:
        out.append(
            (
                stain,
                bids(
                    root=root,
                    datatype="tabular",
                    seg="{seg}",
                    from_="{template}",
                    stain=stain,
                    level=config["fibre_output_level"],
                    desc="frangi",
                    suffix="fibrestats.tsv",
                    **inputs["spim"].wildcards,
                ).format(**wildcards),
            )
        )
    return out


def get_coloc_tsv_input(wildcards):
    """Colocalization table, when two or more stains have objects."""
    if not do_coloc:
        return []
    return rules.merge_into_colocsegstats_tsv.output.tsv.format(**wildcards)


def get_coverage_tsv_input(wildcards):
    """Field-of-view coverage of each region, for additional scans (--aux_spim) only."""
    if not is_aux_scan(wildcards):
        return []
    return rules.aux_coverage_tsv.output.tsv.format(**wildcards)


rule merge_indiv_and_coloc_segstats_tsv:
    """One per-region table of every stain, from its primary methods."""
    input:
        coloc_tsv=get_coloc_tsv_input,
        indiv_tsvs=lambda wildcards: [p for _, p in merged_table_sources(wildcards)],
        coverage_tsv=get_coverage_tsv_input,
    output:
        merged_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            suffix="mergedsegstats.tsv",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_coloc"
    threads: 1
    resources:
        mem_mb=1500,
        runtime=10,
    params:
        # aligned with indiv_tsvs
        stains=lambda wildcards: [s for s, _ in merged_table_sources(wildcards)],
    script:
        "../scripts/merge_indiv_and_coloc_segstats_tsv.py"
