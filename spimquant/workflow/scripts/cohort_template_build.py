"""Run FireANTs build_template, also saving each subject's final-iteration warps.

build_template only saves the normalized moved images, so warps are captured by
patching two functions in its registration_pipeline module: register_batch, to
know the current subject and iteration, and add_shape, which is called with the
finished registration object. Forward and inverse warps are written in ANTs
format to {save_dir}/warps/, named after the subject's first image in the image
list.

Takes the same command line as python -m fireants.scripts.template.build_template.
"""

import os
import runpy

from fireants.scripts.template import registration_pipeline

_register_batch = registration_pipeline.register_batch
_add_shape = registration_pipeline.add_shape

# (warp, inverse warp) filenames for the current batch, in the last iteration
_warp_files = None


def register_batch(
    init_template_batch,
    moving_images_batch,
    args,
    logger,
    avg_warp,
    identifiers,
    is_last_epoch=False,
):
    global _warp_files
    # fail in the first iteration, rather than when saving in the last one
    if not args.do_deform:
        raise ValueError("Saving warps requires do_deform")
    if len(identifiers) != 1:
        raise ValueError("Saving warps requires batch_size: 1")

    _warp_files = None
    if is_last_epoch:
        warp_dir = os.path.join(args.save_dir, "warps")
        os.makedirs(warp_dir, exist_ok=True)
        prefix = os.path.join(warp_dir, identifiers[0].removesuffix(".nii.gz"))
        _warp_files = (f"{prefix}_warp.nii.gz", f"{prefix}_invwarp.nii.gz")
    return _register_batch(
        init_template_batch,
        moving_images_batch,
        args,
        logger,
        avg_warp,
        identifiers,
        is_last_epoch,
    )


def add_shape(avg_warp, reg):
    if _warp_files is not None:
        reg.save_as_ants_transforms([_warp_files[0]])
        reg.save_as_ants_transforms([_warp_files[1]], save_inverse=True)
    return _add_shape(avg_warp, reg)


registration_pipeline.register_batch = register_batch
registration_pipeline.add_shape = add_shape

# as __main__, like python -m: hydra then finds the configs relative to the file
runpy.run_module(
    "fireants.scripts.template.build_template", run_name="__main__", alter_sys=True
)
