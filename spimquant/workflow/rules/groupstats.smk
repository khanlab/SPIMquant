"""
Group-level statistical analysis rules for SPIMquant.

This module performs group-based statistical tests on segmentation statistics
(e.g., fieldfrac, density, volume) across participants, using metadata from
participants.tsv to define contrasts.

It also provides the ``group_otsu`` rule which aggregates per-subject intensity
histograms (produced by ``compute_subject_histogram``) to compute a single set
of Otsu thresholds shared across all subjects.
"""


rule group_otsu:
    """Compute group-level Otsu thresholds from aggregated per-subject histograms.

    Collects the intensity histogram NPZ files computed by
    ``compute_subject_histogram`` for all subjects, merges them onto a
    common intensity grid, and applies multi-level Otsu thresholding to the
    aggregate histogram.  The resulting thresholds are saved as a JSON file
    (consumed by ``multiotsu_group`` during participant-level segmentation)
    and as a PNG figure for visual inspection.

    This rule is the target of ``all_group_otsu`` and should be run before
    participant-level segmentation when ``groupotsu+k{}i{}`` is used as the
    segmentation method.
    """
    input:
        histogram_npz=lambda wildcards: inputs["spim"].expand(
            bids(
                root=work,
                datatype="seg",
                stain=wildcards.stain,
                level=wildcards.level,
                desc="groupotsu+k{k}i{i}".format(k=wildcards.k, i=wildcards.i),
                suffix="histogram.npz",
                **inputs["spim"].wildcards,
            )
        ),
    params:
        otsu_k=lambda wildcards: int(wildcards.k),
        otsu_threshold_index=lambda wildcards: int(wildcards.i),
    output:
        thresholds_json=bids(
            root=root,
            datatype="group",
            stain="{stain}",
            level="{level}",
            desc="groupotsu+k{k,[0-9]+}i{i,[0-9]+}",
            suffix="thresholds.json",
        ),
        thresholds_png=bids(
            root=root,
            datatype="group",
            stain="{stain}",
            level="{level}",
            desc="groupotsu+k{k,[0-9]+}i{i,[0-9]+}",
            suffix="thresholds.png",
        ),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=10,
    script:
        "../scripts/group_otsu.py"


rule perform_group_stats:
    """Perform group-based statistical tests on segmentation statistics.
    
    This rule reads segstats.tsv files from all participants and performs
    statistical tests based on contrasts defined in participants.tsv.
    """
    input:
        segstats_tsvs=lambda wildcards: inputs["spim"].expand(
            bids(
                root=root,
                datatype="tabular",
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
        mem_mb=1500,
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
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
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
        runtime=15,
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
        dseg=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
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
        runtime=15,
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
                datatype="tabular",
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
        mem_mb=1500,
        runtime=10,
    script:
        "../scripts/concat_subj_parquet.py"


rule group_counts_per_voxel:
    """Calculate counts per voxel based on concatenated points
    in template space"""
    input:
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
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
            suffix="{stain}+count.nii.gz",
        ),
    threads: 16
    resources:
        mem_mb=200000,
        runtime=30,
    script:
        "../scripts/counts_per_voxel_template.py"


rule group_coloc_counts_per_voxel:
    """Calculate counts per voxel based on concatenated coloc points
    in template space"""
    input:
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
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
            suffix="coloccount.nii.gz",
        ),
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
                datatype="tabular",
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
        mem_mb=1500,
        runtime=10,
    script:
        "../scripts/concat_subj_parquet_contrast.py"


rule group_counts_per_voxel_contrast:
    """Calculate counts per voxel based on concatenated points
    in template space, filtered by contrast"""
    input:
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
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
            suffix="{stain}+count.nii.gz",
        ),
    threads: 16
    resources:
        mem_mb=200000,
        runtime=30,
    script:
        "../scripts/counts_per_voxel_template.py"


rule group_coloc_counts_per_voxel_contrast:
    """Calculate counts per voxel based on concatenated coloc points
    in template space, filtered by contrast"""
    input:
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
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
            suffix="coloccount.nii.gz",
        ),
    threads: 16
    resources:
        mem_mb=15000,
        runtime=10,
    script:
        "../scripts/coloc_per_voxel_template.py"


rule concat_subj_segstats_contrast:
    """Concatenate segstats.tsv files across subjects filtered by contrast
    and compute group averages.
    
    This rule collects mergedsegstats.tsv files from all participants, 
    adds a participant_id column to identify each subject's data, 
    merges with participant metadata from participants.tsv, filters 
    to include only rows where the contrast_column matches the 
    contrast_value, and computes group averages for each atlas region.
    """
    input:
        segstats_tsvs=lambda wildcards: inputs["spim"].expand(
            bids(
                root=root,
                datatype="tabular",
                seg=wildcards.seg,
                from_=wildcards.template,
                desc=wildcards.desc,
                suffix="mergedsegstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
        participants_tsv=os.path.join(config["bids_dir"], "participants.tsv"),
    params:
        contrast_column="{contrast_column}",
        contrast_value="{contrast_value}",
        metric_columns=expand(
            "{stain}+{metric}", stain=stains_for_seg, metric=config["seg_metrics"]
        ),
        coloc_metric_columns=expand(
            "coloc+{metric}", metric=config["coloc_seg_metrics"]
        ),
    output:
        tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            suffix="groupavgsegstats.tsv",
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=10,
    script:
        "../scripts/concat_subj_segstats_contrast.py"


rule map_groupavg_segstats_to_template_nii:
    """Map group-averaged segstats to template space as NIfTI files.
    
    This rule takes the group-averaged segstats for a specific contrast
    and paints brain regions with the averaged metric values to create 
    volumetric maps for 3D visualization.
    """
    input:
        tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            suffix="groupavgsegstats.tsv",
        ),
        dseg=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    params:
        label_column="index",
        feature_column="{metric}",
    output:
        nii=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            space="{template}",
            desc="{desc}",
            contrast="{contrast_column}+{contrast_value}",
            metric="{metric}",
            suffix="groupavg.nii.gz",
        ),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/map_tsv_dseg_to_nii.py"
