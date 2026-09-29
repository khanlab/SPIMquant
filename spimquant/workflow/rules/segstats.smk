
rule map_regionprops_to_atlas_rois:
    input:
        regionprops_parquet=get_regionprops_parquet,
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
    params:
        coord_column_names=config["coord_column_names"],
    output:
        regionprops_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
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
                suffix="countstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=15,
    script:
        "../scripts/map_atlas_to_regionprops.py"


rule map_coloc_to_atlas_rois:
    input:
        coloc_parquet=bids(
            root=root,
            datatype="tabular",
            space="{template}",
            desc="{desc}",
            suffix="coloc.parquet",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    params:
        coord_column_names=["template_coloc_x", "template_coloc_y", "template_coloc_z"],
    output:
        coloc_tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
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
                desc="{desc}",
                suffix="coloccountstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=15,
    script:
        "../scripts/map_atlas_to_coloc.py"


rule merge_into_segstats_tsv:
    input:
        regionprops_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="regionpropstats.tsv",
            **inputs["spim"].wildcards,
        ),
        counts_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="countstats.tsv",
            **inputs["spim"].wildcards,
        ),
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
        tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                stain="{stain}",
                level="{level}",
                desc="{desc}",
                suffix="segstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=15,
    script:
        "../scripts/merge_into_segstats_tsv.py"


def get_coloc_reference_fieldfrac(wildcards):
    """Fieldfrac table used only to recover ROI volume, so count becomes density.

    Colocalization is cross-stain, so every stain this method segmented carries the
    same ROI volume column and the first one is an arbitrary but stable choice.
    Which stains a method segments depends on the method -- the plaque method covers
    only its own stain -- so this cannot be resolved at parse time.
    """
    stains = stains_for_desc(wildcards.desc)
    if not stains:
        raise ValueError(
            f"--seg_method {wildcards.desc} segments no stains, so colocalization "
            f"statistics cannot be computed for it"
        )
    return bids(
        root=root,
        datatype="tabular",
        seg="{seg}",
        from_="{template}",
        stain=stains[0],
        level=config["registration_level"],
        desc="{desc}",
        suffix="fieldfracstats.tsv",
        **inputs["spim"].wildcards,
    ).format(**wildcards)


rule merge_into_colocsegstats_tsv:
    """ also includes fieldfracstats.tsv to obtain the volume to turn count into density"""
    input:
        coloc_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="colocstats.tsv",
            **inputs["spim"].wildcards,
        ),
        counts_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="coloccountstats.tsv",
            **inputs["spim"].wildcards,
        ),
        fieldfrac_tsv=get_coloc_reference_fieldfrac,
    params:
        columns_to_drop=["fieldfrac"],
    output:
        tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                suffix="colocsegstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    script:
        "../scripts/merge_into_segstats_tsv.py"


def get_coloc_tsv_input(wildcards):
    """Colocalization table, only for methods that segment two or more stains.

    Single-stain methods (the plaque ensemble) have nothing to colocalize, so this
    resolves to an empty list and the merge falls back to the per-stain tables.
    """
    if wildcards.desc not in coloc_seg_methods:
        return []
    return bids(
        root=root,
        datatype="tabular",
        seg="{seg}",
        from_="{template}",
        desc="{desc}",
        suffix="colocsegstats.tsv",
        **inputs["spim"].wildcards,
    ).format(**wildcards)


def get_indiv_segstats_tsvs(wildcards):
    """Per-stain segstats tables for the stains this method actually segmented."""
    paths = expand(
        bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            level=config["registration_level"],
            desc="{desc}",
            suffix="segstats.tsv",
            **inputs["spim"].wildcards,
        ),
        stain=stains_for_desc(wildcards.desc),
        allow_missing=True,
    )
    # Input functions are handed to snakemake verbatim, so the remaining wildcards
    # have to be substituted here rather than left for the usual resolution pass.
    return [path.format(**wildcards) for path in paths]


rule merge_indiv_and_coloc_segstats_tsv:
    input:
        coloc_tsv=get_coloc_tsv_input,
        indiv_tsvs=get_indiv_segstats_tsvs,
    params:
        # must stay aligned with indiv_tsvs, which is per-method
        stains=lambda wildcards: stains_for_desc(wildcards.desc),
    output:
        merged_tsv=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="mergedsegstats.tsv",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    script:
        "../scripts/merge_indiv_and_coloc_segstats_tsv.py"
