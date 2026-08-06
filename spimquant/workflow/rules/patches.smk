"""
Patch extraction workflow for SPIMquant.

This module extracts fixed-size 3D patches from SPIM data and segmentation masks
based on atlas regions. Patches are useful for:
- Training machine learning models
- Visual quality control of segmentation
- Creating datasets for downstream analysis
- Extracting high-resolution regions of interest

Key features:
1. Atlas-guided sampling (extract patches from specific brain regions)
2. Supports raw SPIM data, corrected data, and segmentation masks
3. Fixed patch size (configurable, default 256^3)
4. Random sampling with reproducible seeds
5. Output as collections of NIfTI files
6. Imaris dataset export for high-resolution visualization

Patches are named with atlas abbreviation and patch number for easy identification.
"""


rule create_spim_patches:
    """Create patches from SPIM zarr data based on atlas regions.

    This rule extracts fixed-size patches from the SPIM zarr data at locations
    sampled from specified atlas regions. Patches are saved as NIfTI files
    named with the atlas, label abbreviation, and patch number.

    Level in the output file is the resolution of the patches, (e.g. level 0)
    but level in the dseg input is the (downsampled) registration_level.
    """
    input:
        spim=inputs["spim"].path,
        dseg=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    params:
        patch_size=config.get("patch_size", [256, 256, 256]),
        n_patches=config.get("n_patches_per_label", 10),
        patch_labels=config.get("patch_labels", None),
        hires_level=0,  #input is the raw data
        seed=config.get("patch_seed", 42),
        zarrnii_kwargs=zarrnii_in_kwargs,
        patch_uint8=not config.get("no_patch_uint8", False),
    output:
        patches_dir=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                seg="{seg}",
                from_="{template}",
                level="{level}",
                desc="raw",
                suffix="SPIM.patches",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 32
    resources:
        mem_mb=32000,
        runtime=30,
    script:
        "../scripts/create_patches.py"


rule create_mask_patches:
    """Create patches from segmentation mask zarr data based on atlas regions.

    This rule extracts fixed-size patches from the cleaned segmentation mask
    zarr data at locations sampled from specified atlas regions. Patches are
    saved as NIfTI files named with the atlas, label abbreviation, and patch number.
    """
    input:
        mask=bids_oz_in(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    params:
        patch_size=config.get("patch_size", [256, 256, 256]),
        n_patches=config.get("n_patches_per_label", 10),
        patch_labels=config.get("patch_labels", None),
        seed=config.get("patch_seed", 42),
        hires_level=config["segmentation_level"],
        patch_uint8=not config.get("no_patch_uint8", False),
    output:
        patches_dir=directory(
            bids(
                root=root,
                datatype="seg",
                stain="{stain}",
                seg="{seg}",
                from_="{template}",
                level="{level}",
                desc="{desc}",
                suffix="mask.patches",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 32
    resources:
        mem_mb=32000,
        runtime=30,
    script:
        "../scripts/create_patches.py"


rule create_corrected_spim_patches:
    """Create patches from corrected SPIM zarr data based on atlas regions.

    This rule extracts fixed-size patches from the intensity-corrected SPIM
    zarr data at locations sampled from specified atlas regions. Patches are
    saved as NIfTI files named with the atlas, label abbreviation, and patch number.
    """
    input:
        corrected=bids_oz_in(
            root=work,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="corrected{corrmethod}",
            suffix="SPIM.{ext}",
            **inputs["spim"].wildcards,
        ),
        dseg=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    params:
        patch_size=config.get("patch_size", [256, 256, 256]),
        n_patches=config.get("n_patches_per_label", 10),
        patch_labels=config.get("patch_labels", None),
        seed=config.get("patch_seed", 42),
        hires_level=config["segmentation_level"],
        patch_uint8=not config.get("no_patch_uint8", False),
    output:
        patches_dir=directory(
            bids(
                root=root,
                datatype="seg",
                stain="{stain}",
                seg="{seg}",
                from_="{template}",
                level="{level}",
                desc="corrected{corrmethod}",
                suffix="SPIM.patches",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 32
    resources:
        mem_mb=32000,
        runtime=30,
    script:
        "../scripts/create_patches.py"


rule create_imaris_crops:
    """Create high-resolution Imaris datasets from SPIM data based on atlas region bounding boxes.

    This rule extracts crops from SPIM zarr data based on bounding boxes of
    specified atlas regions. Crops are saved as Imaris datasets using zarrnii's
    to_imaris() function. Level defaults to 0 for high-resolution output.
    """
    input:
        spim=inputs["spim"].path,
        dseg=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    params:
        crop_labels=config.get("crop_labels", None),
        hires_level=0,  # input is the raw data
        zarrnii_kwargs=zarrnii_in_kwargs,
    output:
        crops_dir=directory(
            bids(
                root=root,
                datatype="seg",
                seg="{seg}",
                from_="{template}",
                level="{level}",
                desc="crop",
                suffix="SPIM.imaris",
                **inputs["spim"].wildcards,
            )
        ),
    threads: 32
    resources:
        mem_mb=32000,
        runtime=60,
    script:
        "../scripts/create_imaris_crops.py"


rule extract_bbox_crop:
    """Extract a user-defined bounding box crop from SPIM zarr data in template space.

    Reads a bounding box specification from the bbox_config YAML (provided via
    --bbox_config) and resamples the SPIM zarr data at the requested input level
    to the bounding box in template space at the desired isotropic resolution.

    The reference volume in template space is constructed from the bounding box
    centre coordinates (RAS mm) and the target voxel size, then the floating SPIM
    zarr is pull-resampled to that reference using the composite inverse warp
    (template → subject) via zarrnii's block-wise interpolation.

    A special case (use_brain_mask: true) uses the template brain mask to define
    the bounding box extent, enabling whole-brain resampling at arbitrary resolution.
    """

    input:
        spim=inputs["spim"].path,
        xfm_composite_inv=bids(
            root=root,
            datatype="xfm",
            from_="{template}",
            to="subject",
            suffix="xfm.nii.gz",
            **inputs["spim"].wildcards,
        ),
        brain_mask=lambda wildcards: (
            bids(
                root=root,
                template=wildcards.template,
                desc="brain",
                suffix="mask.nii.gz",
            )
            if bbox_cfgs_by_name.get(wildcards.desc, {}).get("use_brain_mask", False)
            else []
        ),
    params:
        bbox_config=lambda wildcards: bbox_cfgs_by_name[wildcards.desc],
        zarrnii_kwargs=zarrnii_in_kwargs,
    output:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            space="{template}",
            res="{res}um",
            desc="{desc}",
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    threads: 32
    resources:
        mem_mb=32000,
        runtime=60,
    script:
        "../scripts/extract_bbox_in_template_space.py"


rule all_bbox_template_crops:
    """Target rule to generate all custom bounding box crops in template space.

    Expands over all bounding boxes defined in the --bbox_config YAML and all
    stains present in the dataset.  Each bounding box produces one NIfTI file
    per subject per stain at the resolution and input level specified in the
    config.
    """

    input:
        [
            inputs["spim"].expand(
                bids(
                    root=root,
                    datatype="micr",
                    stain="{stain}",
                    space="{template}",
                    res=f"{bbox['resolution_um']}um",
                    desc=bbox["name"],
                    level=bbox.get("input_level", 0),
                    suffix="SPIM.nii.gz",
                    **inputs["spim"].wildcards,
                ),
                stain=stains,
                template=config["template"],
            )
            for bbox in bbox_configs
        ],
