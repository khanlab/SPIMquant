"""
Group-level statistical analysis rules for SPIMquant.

This module performs group-based statistical tests on segmentation statistics
(e.g., fieldfrac, density, volume) across participants, using metadata from
participants.tsv to define contrasts.
"""


rule perform_group_stats:
    """Perform group-based statistical tests on segmentation statistics.
    
    This rule reads segstats.tsv files from all participants and performs
    statistical tests based on contrasts defined in participants.tsv.
    """
    input:
        segstats_tsvs=lambda wildcards: expand(
            bids(
                root=root,
                datatype="micr",
                seg=wildcards.seg,
                from_=wildcards.template,
                stain=wildcards.stain,
                level=config["registration_level"],
                desc=wildcards.desc,
                suffix="segstats.tsv",
                **inputs["spim"].wildcards,
            ),
            zip,
            **inputs["spim"].zip_lists,
        ),
        participants_tsv=os.path.join(config["bids_dir"], "participants.tsv"),
    params:
        contrast_column=config.get("contrast_column", None),
        contrast_values=config.get("contrast_values", None),
    output:
        stats_tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="groupstats.tsv",
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=10,
    script:
        "../scripts/perform_group_stats.py"


rule create_stats_heatmap:
    """Create heatmap visualizations from group statistics results.
    
    This rule takes the group statistics TSV and creates heatmaps for
    visualization of significant differences across brain regions.
    """
    input:
        stats_tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="groupstats.tsv",
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        metric_columns=config.get("stats_metric_columns", ["fieldfrac", "density", "count", "volume"]),
    output:
        heatmap_png=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="groupstats.png",
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=5,
    script:
        "../scripts/create_stats_heatmap.py"


rule map_groupstats_to_template_nii:
    """Map group statistics to template space as NIfTI files.
    
    This rule paints brain regions with statistical values (e.g., t-statistics,
    p-values) to create volumetric heatmaps for 3D visualization.
    """
    input:
        stats_tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="groupstats.tsv",
        ),
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        label_column="index",
        feature_column="{metric}",
    output:
        nii=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            space="{template}",
            stain="{stain}",
            desc="{desc}",
            metric="{metric}",
            suffix="groupstats.nii",
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/map_tsv_dseg_to_nii.py"
