
rule map_segstats_tsv_dseg_to_template_nii:
    """ uses generic script that paints regions with column data (e.g. use this to make density heat-maps)"""
    input:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="mergedsegstats.tsv",
            **inputs["spim"].wildcards,
        ),
        dseg=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.nii.gz"
        ),
        label_tsv=bids_tpl(
            root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"
        ),
    params:
        label_column="index",
        feature_column="{suffix}",
    output:
        nii=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            space="{template}",
            desc="{desc}",
            suffix="{suffix}.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    script:
        "../scripts/map_tsv_dseg_to_nii.py"


rule map_segstats_tsv_dseg_to_subject_nii:
    """ uses generic script that paints regions with column data (e.g. use this to make density heat-maps)"""
    input:
        tsv=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="mergedsegstats.tsv",
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
        label_column="index",
        feature_column="{suffix}",
    output:
        nii=bids(
            root=root,
            datatype="micr",
            seg="{seg}",
            level="{level}",
            from_="{template}",
            desc="{desc}",
            suffix="{suffix}.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    script:
        "../scripts/map_tsv_dseg_to_nii.py"


rule deform_fieldfrac_nii_to_template_nii:
    input:
        flo=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="fieldfrac.nii.gz",
            **inputs["spim"].wildcards,
        ),
        ref=rules.import_template_anat.output.anat,
        xfm_itk=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            type_="itk",
            desc="affine",
            suffix="xfm.txt",
            **inputs["spim"].wildcards,
        ),
        warp=bids(
            root=root,
            datatype="warps",
            from_="subject",
            to="{template}",
            suffix="warp.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            space="{template}",
            suffix="fieldfrac.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    resources:
        mem_mb=32000,
        runtime=30,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear "
        " -i {input.flo} -o {output.nii} "
        " -r {input.ref} -t {input.warp} {input.xfm_itk}"
