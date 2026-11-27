# TEMPLATE: Snakemake rules for a model run through the tiled-inference engine.
#
# Copy into a rules file (e.g. spimquant/workflow/rules/<yourmodel>.smk, and
# `include:` it from the Snakefile) and replace `mymodel` throughout.  The
# config keys it reads are listed in docs/howto/deep_learning_models.md.


rule import_mymodel_model:
    """Download the weights once; inference rules read the local copy."""
    input:
        model=storage(config["models"]["mymodel"]),
    output:
        "resources/models/mymodel.pt",
    localrule: True
    shell:
        "cp {input} {output}"


rule run_mymodel:
    """TODO: what the model segments, at what resolution, from which input."""
    input:
        # RAW intensities, unless the model was trained on corrected data.
        spim=spim_input,
        spim_json=spim_json_input,
        model="resources/models/mymodel.pt",
        # Low-res brain mask: blocks and tiles outside it are skipped. Drop
        # this input (and pass mask_path=None) to run on the whole image.
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["correction_level"],
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        mask=bids_oz_out(
            root=root,
            datatype="seg",  # TODO: or a datatype of its own, e.g. "vessels"
            stain="{stain}",
            level="{level}",
            desc="mymodel",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    # TODO: restrict to the stains the model was trained on, so a request for
    # another stain fails when the DAG is built rather than after GPU hours.
    # Per rule, never globally (a module-level constraint applies everywhere).
    wildcard_constraints:
        stain="CD31|Lectin",
    # Blocks in flight; each holds its halo-grown input and merge buffers.
    threads: 8 * config["mymodel_n_gpus"] if gpu_available else 8
    resources:
        gpu=config["mymodel_n_gpus"] if gpu_available else 0,
        # One task that sees every GPU. Without this the SLURM executor asks
        # for one task per GPU, and srun launches the whole job N times.
        tasks_per_gpu=0,
        cpus_per_gpu=8 if gpu_available else 0,
        # TODO: measure and refit. Starting point: ~3 GB per block in flight
        # (512^3 chunk + 64 halo) plus ~1.6 GB per reader process.
        mem_mb=scaled_mem_mb(24000 + 1700 * config["inference_read_procs"], 20),
        # TODO: measure and refit. Inference splits across GPUs; the pyramid
        # write at the end does not.
        runtime=scaled_runtime(30, 0.15 + 1.0 / config["mymodel_n_gpus"]),
    params:
        zarrnii_kwargs=zarrnii_in_kwargs,
        tile=config["mymodel_tile"],
        stride=config["mymodel_stride"],
        batch_size=config["mymodel_batch_size"],
        chunk=config["mymodel_chunk"],
        merge_mode=config["mymodel_merge_mode"],
        threshold=config["mymodel_threshold"],
        n_gpus=config["mymodel_n_gpus"],
        read_procs=config["inference_read_procs"],
    script:
        "../scripts/mymodel_seg.py"
