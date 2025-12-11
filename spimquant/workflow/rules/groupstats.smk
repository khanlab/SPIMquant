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
        segstats_tsvs=lambda wildcards: inputs["spim"].expand(
            bids(
                root=root,
                datatype="micr",
                seg=wildcards.seg,
                from_=wildcards.template,
                desc=wildcards.desc,
                suffix="mergedsegstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
        participants_tsv=os.path.join(config["bids_dir"], "participants.tsv"),
    params:
        contrast_column=config.get("contrast_column", None),
        contrast_values=config.get("contrast_values", None),
        metric_columns=expand(
            "{stain}+{metric}", stain=stains_for_seg, metric=config["seg_metrics"]
        ),
        coloc_metric_columns=expand(
            "coloc+{metric}", metric=config["coloc_seg_metrics"]
        ),
    output:
        stats_tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
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
            desc="{desc}",
            suffix="groupstats.tsv",
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        metric_columns=expand(
            "{stain}+{metric}", stain=stains_for_seg, metric=config["seg_metrics"]
        ),
        coloc_metric_columns=expand(
            "coloc+{metric}", metric=config["coloc_seg_metrics"]
        ),
    output:
        heatmap_png=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
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
        tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
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
        feature_column="{metric}_{stat}",
    output:
        nii=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            space="{template}",
            desc="{desc}",
            metric="{metric}",
            suffix="{stat}.nii.gz",
        ),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=5,
    script:
        "../scripts/map_tsv_dseg_to_nii.py"


rule concat_subj_parquet:
    """Concatenate parquet files across all subjects.
    
    This rule collects regionprops.parquet or coloc.parquet files
    from all participants, adds a participant_id column to 
    identify each subject's data, and merges with participant 
    metadata from participants.tsv.
    """
    input:
        parquet_files=inputs["spim"].expand(
            bids(
                root=root,
                datatype="micr",
                space="{template}",
                desc="{desc}",
                suffix="{suffix}.parquet",
                **inputs["spim"].wildcards,
            ),
            allow_missing=True,
        ),
        participants_tsv=os.path.join(config["bids_dir"], "participants.tsv"),
    output:
        parquet=bids(
            root=root,
            datatype="group",
            space="{template}",
            desc="{desc}",
            suffix="{suffix,regionprops|coloc}.parquet",
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=10,
    script:
        "../scripts/concat_subj_parquet.py"


rule group_counts_per_voxel:
    """Calculate counts per voxel based on concatenated points
    in template space"""
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        regionprops_parquet=bids(
            root=root,
            datatype="group",
            space="{template}",
            desc="{desc}",
            suffix="regionprops.parquet",
        ),
    params:
        coord_column_names=config["template_coord_column_names"],
    output:
        counts_nii=bids(
            root=root,
            datatype="group",
            space="{template}",
            level="{level}",
            desc="{desc}",
            suffix="{stain}+count.nii",
        ),
    group:
        "subj"
    threads: 16
    resources:
        mem_mb=15000,
        runtime=10,
    script:
        "../scripts/counts_per_voxel_template.py"


rule group_coloc_counts_per_voxel:
    """Calculate counts per voxel based on concatenated coloc points
    in template space"""
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        coloc_parquet=bids(
            root=root,
            datatype="group",
            space="{template}",
            desc="{desc}",
            suffix="coloc.parquet",
        ),
    params:
        coord_column_names=config["template_coloc_coord_column_names"],
    output:
        counts_nii=bids(
            root=root,
            datatype="group",
            space="{template}",
            level="{level}",
            desc="{desc}",
            suffix="coloccount.nii",
        ),
    group:
        "subj"
    threads: 16
    resources:
        mem_mb=15000,
        runtime=10,
    script:
        "../scripts/coloc_per_voxel_template.py"


rule concat_subj_parquet_contrast:
    """Concatenate parquet files across subjects filtered by contrast.
    
    This rule collects regionprops.parquet or coloc.parquet files
    from all participants, adds a participant_id column to 
    identify each subject's data, merges with participant 
    metadata from participants.tsv, and filters to include only
    rows where the contrast_column matches the contrast_value.
    """
    input:
        parquet_files=inputs["spim"].expand(
            bids(
                root=root,
                datatype="micr",
                space="{template}",
                desc="{desc}",
                suffix="{suffix}.parquet",
                **inputs["spim"].wildcards,
            ),
            allow_missing=True,
        ),
        participants_tsv=os.path.join(config["bids_dir"], "participants.tsv"),
    params:
        contrast_column="{contrast_column}",
        contrast_value="{contrast_value}",
    output:
        parquet=bids(
            root=root,
            datatype="group",
            space="{template}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            suffix="{suffix,regionprops|coloc}.parquet",
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=10,
    script:
        "../scripts/concat_subj_parquet_contrast.py"


rule group_counts_per_voxel_contrast:
    """Calculate counts per voxel based on concatenated points
    in template space, filtered by contrast"""
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        regionprops_parquet=bids(
            root=root,
            datatype="group",
            space="{template}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            suffix="regionprops.parquet",
        ),
    params:
        coord_column_names=config["template_coord_column_names"],
    output:
        counts_nii=bids(
            root=root,
            datatype="group",
            space="{template}",
            level="{level}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            suffix="{stain}+count.nii",
        ),
    group:
        "subj"
    threads: 16
    resources:
        mem_mb=15000,
        runtime=10,
    script:
        "../scripts/counts_per_voxel_template.py"


rule group_coloc_counts_per_voxel_contrast:
    """Calculate counts per voxel based on concatenated coloc points
    in template space, filtered by contrast"""
    input:
        template=bids_tpl(root=root, template="{template}", suffix="anat.nii.gz"),
        coloc_parquet=bids(
            root=root,
            datatype="group",
            space="{template}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            suffix="coloc.parquet",
        ),
    params:
        coord_column_names=config["template_coloc_coord_column_names"],
    output:
        counts_nii=bids(
            root=root,
            datatype="group",
            space="{template}",
            level="{level}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            suffix="coloccount.nii",
        ),
    group:
        "subj"
    threads: 16
    resources:
        mem_mb=15000,
        runtime=10,
    script:
        "../scripts/coloc_per_voxel_template.py"
