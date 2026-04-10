"""
Group-level statistical analysis rules for SPIMquant.

This module performs formula-based statistical modelling on segmentation
statistics (e.g., fieldfrac, density, volume) across participants, using
metadata from participants.tsv to fit OLS models and compute pairwise contrasts.
"""


rule perform_group_stats:
    """Perform formula-based group statistical tests on segmentation statistics.

Fits a single global OLS model per region/metric using the user-supplied
formula, then computes pairwise contrast statistics (t-stat, p-value,
Cohen's d) for the specified contrast using the model's covariance matrix.
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
    output:
        stats_tsv=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            contrast="{pairwise_contrast}",
            suffix="groupstats.tsv",
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=10,
    params:
        model=config.get("model", None),
        pairwise_contrast_info=lambda wc: pairwise_contrast_info.get(
            wc.pairwise_contrast, {}
        ),
        within_factors=config.get("within") or [],
        metric_columns=expand(
            "{stain}+{metric}", stain=stains_for_seg, metric=config["seg_metrics"]
        ),
        coloc_metric_columns=expand(
            "coloc+{metric}", metric=config["coloc_seg_metrics"]
        ),
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
            contrast="{pairwise_contrast}",
            suffix="groupstats.tsv",
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    output:
        heatmap_png=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            contrast="{pairwise_contrast}",
            suffix="groupstats.png",
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    params:
        metric_columns=expand(
            "{stain}+{metric}", stain=stains_for_seg, metric=config["seg_metrics"]
        ),
        coloc_metric_columns=expand(
            "coloc+{metric}", metric=config["coloc_seg_metrics"]
        ),
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
            contrast="{pairwise_contrast}",
            suffix="groupstats.tsv",
        ),
        dseg=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    output:
        nii=bids(
            root=root,
            datatype="group",
            seg="{seg}",
            space="{template}",
            desc="{desc}",
            contrast="{pairwise_contrast}",
            metric="{metric}",
            suffix="{stat}.nii.gz",
        ),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        label_column="index",
        feature_column="{metric}_{stat}",
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
    params:
        coord_column_names=config["template_coord_column_names"],
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
    params:
        coord_column_names=config["template_coloc_coord_column_names"],
    script:
        "../scripts/coloc_per_voxel_template.py"
