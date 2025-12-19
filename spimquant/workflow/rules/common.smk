from pathlib import Path
from snakebids import bids


def bids_tpl(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
    return str(Path(bids(root=root, tpl=template)) / bids(tpl=template, **entities))


def resources_path(path):
    """Get path relative to the resources folder"""
    if path.startswith(("http://", "https://")):
        return path
    else:
        return str(Path(workflow.basedir).parent / "resources" / path)


def get_template_path(root, template, template_crop=None):
    """Get template path, optionally cropped based on hemisphere"""
    if template_crop is not None:
        return bids_tpl(
            root=root,
            template=template,
            desc=f"{template_crop}crop",
            suffix="anat.nii.gz",
        )
    else:
        return bids_tpl(root=root, template=template, suffix="anat.nii.gz")


def get_template_for_reg(wildcards):
    """Get the appropriate template file for registration, cropped if specified"""
    if config.get("template_crop") is not None:
        return get_template_path(root, wildcards.template, config["template_crop"])
    else:
        return bids_tpl(root=root, template=wildcards.template, suffix="anat.nii.gz")


def get_stains_all_subjects():

    stain_sets = [
        set(ZarrNii.from_ome_zarr(zarr).list_channels())
        for zarr in inputs["spim"].expand()
    ]
    if all(s == stain_sets[0] for s in stain_sets):
        return list(stain_sets[0])
    else:
        raise ValueError(f"stains across subjects are not consistent, {stain_sets}")
