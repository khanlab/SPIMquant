"""
Reuse of outputs from a previous SPIMquant run (--reuse_derivatives).

Points the workflow at an existing SPIMquant output dataset so that expensive
per-subject products (brain masks, template registration transforms) are read
from it instead of being computed again, e.g. register once, then run several
segmentation configurations against the same masks and transforms.

Each reuse rule copies a file from the prior dataset to the path the computing
rule would have written in the current output dataset, and a ruleorder prefers
it over that computing rule. Downstream rules are unaware of where the file came
from. Because the prior dataset is SPIMquant's own output, its paths are built
from the same bids() patterns as the computing rules (rather than queried with
pybids, which does not know entities like level or type).

When a file is absent from the prior dataset (a subject that was not run, or a
different --registration_level or --template), the reuse rule's input is missing
and Snakemake falls back to the computing rule for that file only.

QC of a reused output would only re-check the prior run, so the registration QC
report is skipped for subjects whose transforms are reused (the brain mask has
no QC output of its own).
"""

reuse_root = config.get("reuse_derivatives")

if reuse_root:
    reuse_root = os.path.abspath(os.path.expandvars(reuse_root))
    if not os.path.isdir(reuse_root):
        raise ValueError(f"--reuse_derivatives: '{reuse_root}' is not a directory")
    if os.path.realpath(reuse_root) == os.path.realpath(root):
        raise ValueError(
            "--reuse_derivatives must be a different dataset from the output directory"
        )
    if not os.path.exists(os.path.join(reuse_root, "dataset_description.json")):
        print(
            f"Warning: --reuse_derivatives '{reuse_root}' has no "
            "dataset_description.json; is it a SPIMquant output dataset?"
        )

reuse_outputs = config.get("reuse_outputs") or []

# --reuse_symlink links to the prior dataset instead of copying (composite warps
# are ~1.5 GB per subject); input paths are absolute, so the links are too
reuse_cmd = "ln -s" if config.get("reuse_symlink") else "cp"


def reuse_path(path):
    """Path in the prior dataset corresponding to an output path under root."""
    return os.path.join(reuse_root, os.path.relpath(str(path), root))


def print_reuse_summary(name, path, **kwargs):
    """Report how many subjects a reuse rule will find files for."""
    paths = inputs["spim"].expand(reuse_path(path), **kwargs)
    n_found = sum(os.path.exists(p) for p in paths)
    print(
        f"Reusing {name} from {reuse_root} for {n_found} of {len(paths)} "
        "subject(s); the rest will be computed"
    )


if reuse_root and "brainmask" in reuse_outputs:

    rule reuse_brainmask:
        """Copy (or link) a brain mask from the prior SPIMquant dataset."""
        input:
            mask=reuse_path(rules.create_mask_from_gmm_and_prior.output.mask),
        output:
            mask=rules.create_mask_from_gmm_and_prior.output.mask,
        localrule: True
        shell:
            f"{reuse_cmd} {{input.mask}} {{output.mask}}"

    ruleorder: reuse_brainmask > create_mask_from_gmm_and_prior > reslice_mask_to_level

    print_reuse_summary(
        "brain masks",
        rules.create_mask_from_gmm_and_prior.output.mask,
        stain=stain_for_reg,
        level=config["registration_level"],
    )


if reuse_root and "transforms" in reuse_outputs:

    rule reuse_composite_warps:
        """Copy (or link) the composite subject<->template warps from the prior SPIMquant dataset."""
        input:
            xfm_composite=reuse_path(
                rules.compose_subject_to_template_warp.output.xfm_composite
            ),
            xfm_composite_inv=reuse_path(
                rules.compose_subject_to_template_warp.output.xfm_composite_inv
            ),
        output:
            xfm_composite=rules.compose_subject_to_template_warp.output.xfm_composite,
            xfm_composite_inv=rules.compose_subject_to_template_warp.output.xfm_composite_inv,
        localrule: True
        shell:
            f"{reuse_cmd} {{input.xfm_composite}} {{output.xfm_composite}} && "
            f"{reuse_cmd} {{input.xfm_composite_inv}} {{output.xfm_composite_inv}}"

    ruleorder: reuse_composite_warps > compose_subject_to_template_warp

    print_reuse_summary(
        "composite warps",
        rules.compose_subject_to_template_warp.output.xfm_composite,
        template=config["template"],
    )


def get_regqc_targets(wildcards):
    """Registration QC reports, for subjects whose transforms are not reused."""
    regqc = inputs["spim"].expand(
        bids(
            root=root,
            datatype="xfm",
            stain="{stain}",
            from_="SPIM",
            to="{template}",
            suffix="regqc.html",
            **inputs["spim"].wildcards,
        ),
        stain=stain_for_reg,
        template=config["template"],
    )
    if not (reuse_root and "transforms" in reuse_outputs):
        return regqc
    # same entity order as regqc, so the lists line up per subject
    reused = inputs["spim"].expand(
        reuse_path(rules.compose_subject_to_template_warp.output.xfm_composite),
        template=config["template"],
    )
    return [qc for qc, xfm in zip(regqc, reused) if not os.path.exists(xfm)]
