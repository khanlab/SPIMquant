rule import_lantern_abeta_fold:
    """Download one fold of the LANTERN Abeta plaque ensemble."""
    input:
        model=lambda wildcards: storage(
            config["models"]["lantern_abeta"][f"fold{wildcards.fold}"]
        ),
    output:
        "resources/models/lantern-ki3-abeta/fold{fold}/seg_model.pt",
    wildcard_constraints:
        fold="[0-9]+",
    localrule: True
    shell:
        "cp {input} {output}"


rule run_lantern_plaques:
    """Segment Abeta plaques with the 5-fold LANTERN ensemble.

    Runs at config['plaque_level'] on RAW intensities -- the model was not
    trained on N4-corrected data, so this bypasses the bias-field chain that
    the GMM and Otsu methods depend on.

    Tiles at 128^3 with stride 64 and combines overlapping tiles using the
    configured merge mode (`max` by default, to match the published LANTERN
    inference). The low-res brain mask restricts inference to blocks that touch
    tissue, which is what makes whole-brain 5-fold inference tractable.

    Writes votes/n_folds rather than a binary mask, so the vote threshold can be
    changed downstream without re-running inference.
    """
    input:
        spim=spim_input,
        models=expand(
            "resources/models/lantern-ki3-abeta/fold{fold}/seg_model.pt",
            fold=range(config["plaque_n_folds"]),
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["correction_level"],
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
    params:
        zarrnii_kwargs=zarrnii_in_kwargs,
        tile=config["plaque_tile"],
        stride=config["plaque_stride"],
        merge_mode=config.get("plaque_merge_mode", "max"),
        batch_size=config["plaque_batch_size"],
        chunk=config["plaque_chunk"],
        n_gpus=config["plaque_n_gpus"],
    output:
        probseg=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc=config["plaque_seg_method"],
            suffix="probseg.{ext}",
            **inputs["spim"].wildcards,
        ),
    wildcard_constraints:
        # Trained on amyloid-beta only. Constraining the stain per-rule (never
        # globally -- a module-level constraint in an included smk applies to the
        # whole workflow) means a target asking for a lantern mask on any other
        # stain fails at DAG build instead of quietly burning a GPU on it.
        stain=stain_for_plaques or "^$",
    threads: 8 * config["plaque_n_gpus"]
    resources:
        gpu=config["plaque_n_gpus"],
        # Suppress the executor's default --ntasks-per-gpu=1. With --gpus=N that
        # asks SLURM for N tasks, so srun launches the whole job N times, each task
        # pinned to one GPU: the ensemble runs N times over the same volume, every
        # copy sees a single device, and they race to write the same output store.
        # Setting this to 0 skips the flag entirely (submit_string.py:82), leaving
        # one task that sees all N GPUs -- which is what DevicePool expects.
        tasks_per_gpu=0,
        cpus_per_gpu=8,
        mem_mb=256000,
        disk_mb=2097152,
        runtime=lambda wildcards: max(
            1,
            int(
                2880.0
                / (3.0 ** float(wildcards.level))
                / float(config["plaque_n_gpus"])
            ),
        ),
        # stride-64 tiling x 5 folds, split across GPUs
    script:
        "../scripts/lantern_plaques.py"


rule binarize_lantern_plaques:
    """Threshold the LANTERN vote map and upsample it to the segmentation level.

    Emits 0/100 at config['segmentation_level'], matching the segmentation.smk
    mask convention, so the standard fieldfrac / regionprops / counts / segstats
    chain consumes it unchanged.
    """
    input:
        probseg=bids_oz_in(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=max(config["plaque_level"],config["segmentation_level"]),
            desc=config["plaque_seg_method"],
            suffix="probseg.{ext}",
            **inputs["spim"].wildcards,
        ),
        spim=spim_input,
    params:
        zarrnii_kwargs=zarrnii_in_kwargs,
        n_folds=config["plaque_n_folds"],
        vote_threshold=config["plaque_vote_threshold"],
        target_level=config["segmentation_level"],
    output:
        mask=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc=config["plaque_seg_method"],
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    wildcard_constraints:
        stain=stain_for_plaques or "^$",
    threads: 32
    resources:
        mem_mb=64000,
        disk_mb=2097152,
        runtime=180,
    script:
        "../scripts/binarize_lantern_plaques.py"
