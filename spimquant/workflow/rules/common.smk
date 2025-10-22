from pathlib import Path
from snakebids import bids


def bids_tpl(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
    return str(Path(bids(root=root, tpl=template)) / bids(tpl=template, **entities))

def resources_path(path):
    """ Get path relative to the resources folder """
    return str(Path(workflow.basedir).parent / 'resources' / path)


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
