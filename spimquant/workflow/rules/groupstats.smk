"""
Group-level statistical analysis rules for SPIMquant.

This module performs formula-based statistical modelling on segmentation
statistics (e.g., fieldfrac, density, volume) across participants, using
metadata from participants.tsv to fit OLS models and compute pairwise contrasts.

Group-level outputs take the reference scans only: additional scans (--aux_spim)
are of the same participants, so would be counted twice. The additional scans
get their own group-level outputs instead, for each acquisition (see
aux_group_entities), with regions their field of view covers too little of left
out of the statistics (--aux_min_coverage).
"""


def get_bids_root_file(filename):
    """Return a file path located at the root of the input BIDS dataset."""
    return os.path.join(str(config["bids_dir"]), filename)


rule perform_group_stats:
    """Perform formula-based group statistical tests on segmentation statistics.

    Fits a single global OLS model per region/metric using the user-supplied
    formula, then computes pairwise contrast statistics (t-stat, p-value,
    Cohen's d) for the specified contrast using the model's covariance matrix.
    """
    input:
        segstats_tsvs=lambda wildcards: ref_spim.expand(
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
        participants_tsv=get_bids_root_file("participants.tsv"),
    output:
        stats_tsv=bids(
            root=group_stats_root,
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            contrast="{pairwise_contrast}",
            suffix="groupstats.tsv",
        ),
    log:
        bids(
            root="logs",
            datatype="group_stats",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            contrast="{pairwise_contrast}",
            suffix="groupstats.log",
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=10,
    params:
        model=config.get("group_stats_model", None),
        pairwise_contrast_info=lambda wc: pairwise_contrast_info.get(
            wc.pairwise_contrast, {}
        ),
        within_factors=config.get("group_stats_within") or [],
        group_stats_where=config.get("group_stats_where", None),
        # only applies to tables with a coverage column (additional scans)
        min_coverage=config["aux_min_coverage"],
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
            root=group_stats_root,
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            contrast="{pairwise_contrast}",
            suffix="groupstats.tsv",
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    output:
        heatmap_png=bids(
            root=group_stats_root,
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
            root=group_stats_root,
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
            root=group_stats_root,
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
        parquet_files=ref_spim.expand(
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
        participants_tsv=get_bids_root_file("participants.tsv"),
    output:
        parquet=bids(
            root=group_stats_root,
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
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
        regionprops_parquet=bids(
            root=group_stats_root,
            space="{template}",
            desc="{desc}",
            suffix="regionprops.parquet",
        ),
    output:
        counts_nii=bids(
            root=group_stats_root,
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
            root=group_stats_root,
            space="{template}",
            desc="{desc}",
            suffix="coloc.parquet",
        ),
    output:
        counts_nii=bids(
            root=group_stats_root,
            space="{template}",
            level="{level}",
            desc="{desc}",
            suffix="coloccount.nii.gz",
        ),
    threads: 16
    resources:
        mem_mb=16000,
        runtime=10,
    params:
        coord_column_names=config["template_coloc_coord_column_names"],
    script:
        "../scripts/coloc_per_voxel_template.py"


rule concat_subj_segstats:
    """Concatenate subject-level segstats TSV files across all participants.

    Merges all individual per-subject mergedsegstats TSV files into a single
    group-level TSV, adding a participant_id column and joining with participant
    metadata from participants.tsv. This file is always produced at group
    analysis level, regardless of whether contrasts are requested, so users can
    export it to their own statistics tools.
    """
    input:
        tsv_files=lambda wildcards: ref_spim.expand(
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
        participants_tsv=get_bids_root_file("participants.tsv"),
    output:
        merged_tsv=bids(
            root=group_stats_root,
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="allsubjects.tsv",
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=10,
    script:
        "../scripts/concat_subj_segstats.py"


rule concat_subj_regionpropstats:
    """Concatenate subject-level regionpropstats TSV files across participants.

    Merges all individual per-subject per-stain regionpropstats TSV files into a
    single group-level TSV, adding a participant_id column and joining with
    participant metadata from participants.tsv. This preserves atlas labels and
    region-level summary statistics for downstream export and QC.
    """
    input:
        tsv_files=lambda wildcards: ref_spim.expand(
            bids(
                root=root,
                datatype="tabular",
                seg=wildcards.seg,
                from_=wildcards.template,
                stain=wildcards.stain,
                level=wildcards.level,
                desc=wildcards.desc,
                suffix="regionpropstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
        participants_tsv=get_bids_root_file("participants.tsv"),
    output:
        merged_tsv=bids(
            root=group_stats_root,
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="allsubjects_regionpropstats.tsv",
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=10,
    script:
        "../scripts/concat_subj_segstats.py"


# Group-level outputs of additional scans (--aux_spim): the rules above, for the
# scans of each acquisition, with its entities (aux_group_entities) in the names
# (only defined with additional scans, whose names would otherwise clash).
if aux_groups:

    def _aux_mergedsegstats(wildcards):
        return aux_group_scans(
            wildcards,
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                suffix="mergedsegstats.tsv",
                **inputs["spim"].wildcards,
            ),
        )

    use rule perform_group_stats as perform_group_stats_aux with:
        input:
            segstats_tsvs=_aux_mergedsegstats,
            participants_tsv=get_bids_root_file("participants.tsv"),
        output:
            stats_tsv=bids(
                root=group_stats_root,
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                contrast="{pairwise_contrast}",
                suffix="groupstats.tsv",
                **aux_group_entities,
            ),
        log:
            bids(
                root="logs",
                datatype="group_stats",
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                contrast="{pairwise_contrast}",
                suffix="groupstats.log",
                **aux_group_entities,
            ),

    use rule create_stats_heatmap as create_stats_heatmap_aux with:
        input:
            stats_tsv=rules.perform_group_stats_aux.output.stats_tsv,
            label_tsv=bids(
                root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
            ),
        output:
            heatmap_png=bids(
                root=group_stats_root,
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                contrast="{pairwise_contrast}",
                suffix="groupstats.png",
                **aux_group_entities,
            ),

    use rule map_groupstats_to_template_nii as map_groupstats_to_template_nii_aux with:
        input:
            tsv=rules.perform_group_stats_aux.output.stats_tsv,
            dseg=bids(
                root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
            ),
            label_tsv=bids(
                root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
            ),
        output:
            nii=bids(
                root=group_stats_root,
                seg="{seg}",
                space="{template}",
                desc="{desc}",
                contrast="{pairwise_contrast}",
                metric="{metric}",
                suffix="{stat}.nii.gz",
                **aux_group_entities,
            ),

    use rule concat_subj_segstats as concat_subj_segstats_aux with:
        input:
            tsv_files=_aux_mergedsegstats,
            participants_tsv=get_bids_root_file("participants.tsv"),
        output:
            merged_tsv=bids(
                root=group_stats_root,
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                suffix="allsubjects.tsv",
                **aux_group_entities,
            ),
