
rule map_regionprops_to_atlas_rois:
    input:
        regionprops_parquet=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            desc="deform",
            level="{level}",
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        coord_column_names=config["coord_column_names"],
    output:
        regionprops_tsv=temp(
            bids(
                root=root,
                datatype="micr",
                seg="{seg}",
                from_="{template}",
                stain="{stain}",
                level="{level}",
                desc="{desc}",
                suffix="regionpropstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
        counts_tsv=temp(
            bids(
                root=root,
                datatype="micr",
                seg="{seg}",
                from_="{template}",
                stain="{stain}",
                level="{level}",
                desc="{desc}",
                suffix="countstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    script:
        "../scripts/map_atlas_to_regionprops.py"


rule map_coloc_to_atlas_rois:
    input:
        coloc_parquet=bids(
            root=root,
            datatype="micr",
            space="{template}",
            desc="{desc}",
            suffix="coloc.parquet",
            **inputs["spim"].wildcards,
        ),
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        coord_column_names=["template_coloc_x", "template_coloc_y", "template_coloc_z"],
    output:
        coloc_tsv=temp(
            bids(
                root=root,
                datatype="micr",
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
                datatype="micr",
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                suffix="coloccountstats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    script:
        "../scripts/map_atlas_to_coloc.py"


rule merge_into_segstats_tsv:
    input:
        regionprops_tsv=bids(
            root=root,
            datatype="micr",
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
            datatype="micr",
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
            datatype="micr",
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
                datatype="micr",
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
        mem_mb=1500,
        runtime=15,
    script:
        "../scripts/merge_into_segstats_tsv.py"


rule merge_into_colocsegstats_tsv:
    """ also includes fieldfracstats.tsv to obtain the volume to turn count into density"""
    input:
        coloc_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="colocstats.tsv",
            **inputs["spim"].wildcards,
        ),
        counts_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="coloccountstats.tsv",
            **inputs["spim"].wildcards,
        ),
        fieldfrac_tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            stain=stains_for_seg[0],
            level=config["registration_level"],
            desc="{desc}",
            suffix="fieldfracstats.tsv",
            **inputs["spim"].wildcards,
        ),
    params:
        columns_to_drop=["fieldfrac"],
    output:
        tsv=temp(
            bids(
                root=root,
                datatype="micr",
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

def get_coloc_tsv_input_kwargs():
    """return coloc_tsv only if we have multiple stains to segment"""
    if len(stains_for_seg) == 1:
        return {}
    else:
        return {
            "coloc_tsv": bids(
                root=root,
                datatype="micr",
                seg="{seg}",
                from_="{template}",
                desc="{desc}",
                suffix="colocsegstats.tsv",
                **inputs["spim"].wildcards,
            )
        }


rule merge_indiv_and_coloc_segstats_tsv:
    input:
        **get_coloc_tsv_input_kwargs(),
        indiv_tsvs=expand(
            bids(
                root=root,
                datatype="micr",
                seg="{seg}",
                from_="{template}",
                stain="{stain}",
                level=config["registration_level"],
                desc="{desc}",
                suffix="segstats.tsv",
                **inputs["spim"].wildcards,
            ),
            stain=stains_for_seg,
            allow_missing=True,
        ),
    params:
        stains=stains_for_seg,
    output:
        merged_tsv=bids(
            root=root,
            datatype="micr",
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


