from pathlib import Path
from snakebids import bids


def bids_tpl(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
    return str(Path(bids(root=root, tpl=template)) / bids(tpl=template, **entities))


def get_template_path(root, template, template_crop=None):
    """Get template path, optionally cropped based on hemisphere"""
    if template_crop is not None:
        return bids_tpl(root=root, template=template, desc=f"{template_crop}crop", suffix="anat.nii.gz")
    else:
        return bids_tpl(root=root, template=template, suffix="anat.nii.gz")
