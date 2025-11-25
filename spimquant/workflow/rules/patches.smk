
rule create_spim_patches:
    """Create patches from SPIM zarr data based on atlas regions.

    This rule extracts fixed-size patches from the SPIM zarr data at locations
    sampled from specified atlas regions. Patches are saved as NIfTI files
    named with the atlas, label abbreviation, and patch number.
    """
    input:
        spim=inputs["spim"].path,
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
        patch_size=config.get("patch_size", [256, 256, 256]),
        n_patches=config.get("n_patches_per_label", 10),
        patch_labels=config.get("patch_labels", None),
        seed=config.get("patch_seed", 42),
        zarrnii_kwargs={"orientation": config["orientation"]},
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
    group:
        "subj"
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
        mask=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            dslevel=config["registration_level"],
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ome.zarr",
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
        patch_size=config.get("patch_size", [256, 256, 256]),
        n_patches=config.get("n_patches_per_label", 10),
        patch_labels=config.get("patch_labels", None),
        seed=config.get("patch_seed", 42),
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        patches_dir=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                seg="{seg}",
                from_="{template}",
                level="{level}",
                desc="{desc}",
                suffix="mask.patches",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj"
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
        corrected=bids(
            root=work,
            datatype="micr",
            stain="{stain}",
            dslevel=config["registration_level"],
            level=config["segmentation_level"],
            desc="corrected",
            corrmethod="{corrmethod}",
            suffix="SPIM.ome.zarr",
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
        patch_size=config.get("patch_size", [256, 256, 256]),
        n_patches=config.get("n_patches_per_label", 10),
        patch_labels=config.get("patch_labels", None),
        seed=config.get("patch_seed", 42),
        zarrnii_kwargs={"orientation": config["orientation"]},
    output:
        patches_dir=directory(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                seg="{seg}",
                from_="{template}",
                level="{level}",
                desc="corrected{corrmethod}",
                suffix="SPIM.patches",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj"
    threads: 32
    resources:
        mem_mb=32000,
        runtime=30,
    script:
        "../scripts/create_patches.py"
