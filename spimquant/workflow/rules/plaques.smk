def lantern_model_spec():
    """(name, revision) from --lantern_config, given as ``<name>[@<revision>]``.

    ``name`` picks the Hugging Face repo ``apooladi/lantern-<name>``; every
    repo under that prefix shares the lantern-ki3-abeta layout, so nothing
    here is specific to one ensemble. ``revision`` defaults to ``main``.
    """
    name, _, revision = str(config["lantern_config"]).partition("@")
    if not name:
        raise ValueError(
            f"--lantern_config {config['lantern_config']!r}: expected <name>[@<revision>]"
        )
    return name, (revision or "main")


lantern_name, lantern_revision = lantern_model_spec()

lantern_fold_path = "resources/models/lantern-{name}/{revision}/fold{fold}/seg_model.pt"


rule import_lantern_fold:
    """Download one fold of a LANTERN plaque ensemble.

    Wildcard-driven so that switching --lantern_config (or its revision) is a
    different target, never a stale cached checkpoint under the same path.
    """
    input:
        model=lambda wildcards: storage(
            config["models"]["lantern"].format(
                name=wildcards.name, revision=wildcards.revision, fold=wildcards.fold
            )
        ),
    output:
        lantern_fold_path,
    wildcard_constraints:
        name="[A-Za-z0-9_.-]+",
        revision="[A-Za-z0-9_.-]+",
        fold="[0-9]+",
    localrule: True
    shell:
        "cp {input} {output}"


rule run_lantern_plaques:
    """Segment Abeta plaques with a 5-fold LANTERN ensemble (--lantern_config).

    Runs on the grid nearest config['plaque_iso_res'] um isotropic (the scale
    the ensemble was trained at), on RAW intensities -- the model was not
    trained on N4-corrected data, so this bypasses the bias-field chain that
    the GMM and Otsu methods depend on. config['plaque_level'] is the starting
    pyramid level; on pyramids that downsample z as well as x/y the script
    loads a finer level and downsamples per axis to land near-isotropic. The
    output is therefore named for the grid (level-neariso4), not a level.

    The probseg has one channel per predicted class beyond background: for a
    binary ensemble just the plaque vote fraction; for a 3-class ensemble a
    second channel with the fraction of folds calling class 2 (the
    "annotated false positive" class -- bright non-plaque structure).

    Tiles at 128^3 with stride 64 and combines overlapping tiles by max over
    votes. The low-res brain mask restricts inference to blocks that touch
    tissue, which is what makes whole-brain 5-fold inference tractable. The
    folds share one frozen encoder, so it runs once per batch and only the
    five decoders are evaluated separately (see lantern_plaques.py).

    Writes votes/n_folds rather than a binary mask, so the vote threshold can be
    changed downstream without re-running inference.
    """
    input:
        spim=spim_input,
        models=expand(
            lantern_fold_path,
            name=lantern_name,
            revision=lantern_revision,
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
        start_level=config["plaque_level"],
        iso_res=config["plaque_iso_res"],
        tile=config["plaque_tile"],
        stride=config["plaque_stride"],
        batch_size=config["plaque_batch_size"],
        chunk=config["plaque_chunk"],
        n_gpus=config["plaque_n_gpus"],
    output:
        probseg=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=plaque_level_label,
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
                / (3.0 ** float(config["plaque_level"]))
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
    chain consumes it unchanged. Only the plaque channel of the probseg reaches
    this mask; the class-2 channel of a 3-class ensemble is used, when
    --label_filter is set, to drop plaque components that touch a large
    class-2 component, and --min_size drops plaque components smaller than
    that many near-iso voxels. Both filters are decided on the near-iso grid
    (connected components there, 26-connectivity) before upsampling.
    """
    input:
        probseg=bids_oz_in(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=plaque_level_label,
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
        min_size=config["min_size"],
        label_filter=config["label_filter"],
        label_filter_size=config["label_filter_size"],
        fp_min_votes=config["plaque_fp_min_votes"],
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
        # Component filtering reads the near-iso probseg once in 512^3 blocks
        # (one per thread in flight, ~1.3 GB each with labels) and keeps only
        # the sparse plaque voxels. Measured on a 45 G-voxel grid with 32
        # threads: 5 min / 42 GB peak for a binary ensemble, 12 min / 63 GB
        # for a 3-class one with --label_filter, before the upsample/write.
        mem_mb=120000,
        disk_mb=2097152,
        runtime=240,
    script:
        "../scripts/binarize_lantern_plaques.py"
