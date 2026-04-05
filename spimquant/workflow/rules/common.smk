"""
Common utility functions for SPIMquant workflows.

This module provides shared helper functions used across multiple workflow rules:

- resources_path(): Path resolver for resources directory
- get_template_path(): Template file locator with optional cropping
- get_template_for_reg(): Registration-specific template selector
- get_stains_all_subjects(): Stain consistency validator across subjects

These functions abstract away path complexity and ensure consistent file
organization following BIDS conventions.
"""

from pathlib import Path


def resources_path(path):
    """Get path relative to the resources folder"""
    if path.startswith(("http://", "https://")):
        return path
    else:
        return str(Path(workflow.basedir).parent / "resources" / path)


def get_template_path(root, template, template_crop=None):
    """Get template path, optionally cropped based on hemisphere"""
    if use_spim_template:
        suffix = stain_for_reg
    else:
        suffix = "anat"

    if template_crop is not None:
        return bids(
            root=root,
            template=template,
            desc=f"{template_crop}crop",
            suffix=f"{suffix}.nii.gz",
        )
    else:
        return bids(root=root, template=template, suffix=f"{suffix}.nii.gz")


def get_template_for_reg(wildcards):
    """Get the appropriate template file for registration, cropped if specified"""
    if use_spim_template:
        suffix = stain_for_reg
    else:
        suffix = "anat"

    if config.get("template_crop") is not None:
        return get_template_path(root, wildcards.template, config["template_crop"])
    else:
        return bids(root=root, template=wildcards.template, suffix=f"{suffix}.nii.gz")


def get_stains_all_subjects():

    stain_sets = [
        set(ZarrNii.from_ome_zarr(zarr).list_channels())
        for zarr in inputs["spim"].expand()
    ]
    if all(s == stain_sets[0] for s in stain_sets):
        return list(stain_sets[0])
    else:
        raise ValueError(f"stains across subjects are not consistent, {stain_sets}")


def get_regionprops_parquet(wildcards):
    if do_vessels:
        return bids(
            root=root,
            datatype="tabular",
            desc="{desc}+vessels",
            space="{template}",
            vessels=stain_for_vessels,
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        )
    else:
        return bids(
            root=root,
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        )


def get_composite_warp(wildcards):
    """Get the composite warp from subject to template space.

    Returns the path to either the direct composite warp (no via-template)
    or the combined via-template composite warp if --via-template is configured.
    """
    spim_wildcards = {k: getattr(wildcards, k) for k in inputs["spim"].wildcards.keys()}
    if via_template is None:
        return bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to=wildcards.template,
            suffix="xfm.nii.gz",
            **spim_wildcards,
        )
    else:
        return bids(
            root=root,
            datatype="xfm",
            from_="subject",
            to=wildcards.template,
            via=via_template,
            suffix="xfm.nii.gz",
            **spim_wildcards,
        )


def get_composite_warp_inv(wildcards):
    """Get the inverse composite warp from template to subject space.

    Returns the path to either the direct inverse composite warp (no via-template)
    or the combined via-template inverse composite warp if --via-template is configured.
    """
    spim_wildcards = {k: getattr(wildcards, k) for k in inputs["spim"].wildcards.keys()}
    if via_template is None:
        return bids(
            root=root,
            datatype="xfm",
            from_=wildcards.template,
            to="subject",
            suffix="xfm.nii.gz",
            **spim_wildcards,
        )
    else:
        return bids(
            root=root,
            datatype="xfm",
            from_=wildcards.template,
            to="subject",
            via=via_template,
            suffix="xfm.nii.gz",
            **spim_wildcards,
        )


def get_affine_warped_for_qc(wildcards):
    """Get the affine-warped SPIM image for QC report.

    Returns the affine-warped image in the appropriate space:
    - Without via-template: image is warped to the final template space
    - With via-template: image is warped to the intermediate (via) template space,
      since that is where the affine registration is performed.
    """
    spim_wildcards = {k: getattr(wildcards, k) for k in inputs["spim"].wildcards.keys()}
    target_space = via_template if via_template is not None else wildcards.template
    return bids(
        root=root,
        datatype="xfm",
        space=target_space,
        stain=wildcards.stain,
        desc="affinewarped",
        suffix="SPIM.nii.gz",
        **spim_wildcards,
    )
