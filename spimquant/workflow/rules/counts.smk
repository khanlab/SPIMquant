rule counts_per_voxel:
    """Calculate counts per voxel based on points"""
    input:
        ref_spim=inputs["spim"].path,
        regionprops_parquet=bids(
            root=root,
            datatype="tabular",
            stain="{stain}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    params:
        coord_column_names=config["coord_column_names"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        counts_nii=bids(
            root=root,
            datatype="seg",
            level="{level}",
            stain="{stain}",
            desc="{desc}",
            suffix="counts.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 16
    resources:
        mem_mb=15000,
        runtime=20,
    script:
        "../scripts/counts_per_voxel.py"


rule counts_per_voxel_template:
    """Calculate counts per voxel based on points
    in template space"""
    input:
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
        regionprops_parquet=bids(
            root=root,
            datatype="tabular",
            space="{template}",
            desc="{desc}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    params:
        coord_column_names=config["template_coord_column_names"],
    output:
        counts_nii=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            space="{template}",
            desc="{desc}",
            suffix="counts.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 16
    resources:
        mem_mb=64000,
        runtime=30,
    script:
        "../scripts/counts_per_voxel_template.py"


rule coloc_per_voxel_template:
    """Calculate coloc counts per voxel based on points
    in template space"""
    input:
        template=bids(root=root, template="{template}", suffix="anat.nii.gz"),
        coloc_parquet=bids(
            root=root,
            datatype="tabular",
            space="{template}",
            desc="{desc}",
            suffix="coloc.parquet",
            **inputs["spim"].wildcards,
        ),
    params:
        coord_column_names=config["template_coloc_coord_column_names"],
    output:
        counts_nii=bids(
            root=root,
            datatype="seg",
            space="{template}",
            desc="{desc}",
            suffix="coloccounts.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 16
    resources:
        mem_mb=64000,
        runtime=30,
    script:
        "../scripts/coloc_per_voxel_template.py"
