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

import json
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


def get_stains_all_subjects(ignore_stains=None):
    """Get stains across all subjects, optionally filtering out ignored stains."""
    ignore_set = set(ignore_stains) if ignore_stains else set()

    stain_sets = []
    for zarr in inputs["spim"].expand():
        channels = set(get_spim_channels(zarr))
        # Remove ignored stains
        channels = channels - ignore_set
        stain_sets.append(channels)

    if all(s == stain_sets[0] for s in stain_sets):
        return sorted(stain_sets[0])
    else:
        raise ValueError(f"stains across subjects are not consistent: {stain_sets}")


def get_spim_json_path(spim_path):
    spim_path = str(spim_path)
    if spim_path.endswith(".ozx"):
        return spim_path[: -len(".ozx")] + ".json"
    if spim_path.endswith(".ome.zarr.zip"):
        return spim_path[: -len(".ome.zarr.zip")] + ".json"
    if spim_path.endswith(".ome.zarr"):
        return spim_path[: -len(".ome.zarr")] + ".json"
    return str(Path(spim_path).with_suffix(".json"))


def _normalize_channel_labels(value):
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        channels = [str(item) for item in value if isinstance(item, (str, int, float))]
        return channels or None
    return None


def _extract_channels_from_json(metadata):
    """Extract channel labels from sidecar keys, using "SampleStaining" key"""
    channels = _normalize_channel_labels(metadata.get("SampleStaining"))
    if channels:
        return channels

    return None


def get_spim_json_overrides(spim_path):
    json_path = get_spim_json_path(spim_path)
    json_file = Path(json_path)
    if not json_file.exists():
        return {}

    with json_file.open() as fp:
        metadata = json.load(fp)

    overrides = {}
    channels = _extract_channels_from_json(metadata)
    if channels:
        overrides["set_channel_labels"] = channels

    orientation = metadata.get("OrientationStringXYZ")
    if orientation:
        overrides["orientation"] = str(orientation)

    return overrides


def get_spim_channels(spim_path):
    overrides = get_spim_json_overrides(spim_path)
    if "set_channel_labels" in overrides:
        return overrides["set_channel_labels"]

    return ZarrNii.from_file(spim_path).list_channels()


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
