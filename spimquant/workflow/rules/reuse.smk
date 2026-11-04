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

import glob as _glob

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


def _reused_warp(rule_output, legacy_output):
    """Input function for a method's composite warp in the prior dataset.

    Datasets from before registration methods were added only have the warp
    without a desc-, which was always computed by the greedy method.
    """

    def get_path(wildcards):
        path = reuse_path(rule_output).format(**wildcards)
        legacy = reuse_path(legacy_output).format(**wildcards)
        if (
            wildcards.reg == "greedy"
            and not os.path.exists(path)
            and os.path.exists(legacy)
            and not os.path.islink(legacy)
        ):
            return legacy
        return path

    return get_path


if reuse_root and "transforms" in reuse_outputs:

    rule reuse_composite_warps:
        """Copy (or link) a method's composite warps from the prior SPIMquant dataset."""
        input:
            xfm_composite=_reused_warp(
                rules.compose_subject_to_template_warp.output.xfm_composite,
                rules.select_subject_to_template_warp.output.xfm_composite,
            ),
            xfm_composite_inv=_reused_warp(
                rules.compose_subject_to_template_warp.output.xfm_composite_inv,
                rules.select_subject_to_template_warp.output.xfm_composite_inv,
            ),
        output:
            xfm_composite=rules.compose_subject_to_template_warp.output.xfm_composite,
            xfm_composite_inv=rules.compose_subject_to_template_warp.output.xfm_composite_inv,
        wildcard_constraints:
            reg=reg_constraint(),
        localrule: True
        shell:
            f"{reuse_cmd} {{input.xfm_composite}} {{output.xfm_composite}} && "
            f"{reuse_cmd} {{input.xfm_composite_inv}} {{output.xfm_composite_inv}}"

    ruleorder: reuse_composite_warps > compose_subject_to_template_warp

    _reuse_method = config["registration_method"]
    _reuse_found = []
    _reuse_legacy = 0
    for desc_path, legacy_path in zip(
        inputs["spim"].expand(
            reuse_path(
                rules.compose_subject_to_template_warp.output.xfm_composite
            ),
            template=config["template"],
            reg=_reuse_method,
        ),
        inputs["spim"].expand(
            reuse_path(
                rules.select_subject_to_template_warp.output.xfm_composite
            ),
            template=config["template"],
        ),
    ):
        _is_legacy = (
            not os.path.exists(desc_path)
            and _reuse_method == "greedy"
            and os.path.exists(legacy_path)
            and not os.path.islink(legacy_path)
        )
        _reuse_found.append(os.path.exists(desc_path) or _is_legacy)
        _reuse_legacy += _is_legacy
    print(
        f"Reusing desc-{_reuse_method} composite warps from {reuse_root} for "
        f"{sum(_reuse_found)} of {len(_reuse_found)} subject(s); the rest will be "
        "computed"
    )
    if _reuse_legacy:
        print(
            f"Warning: {_reuse_legacy} reused warp(s) have no desc- entity, so they "
            "predate registration methods and may have been registered to a packaged "
            "SPIM template, which is not accurately aligned to the template's atlases. "
            "Atlas labels from these warps will inherit that error: use them only as "
            "the starting point for a cohort template, or register again."
        )


def get_regqc_targets(wildcards):
    """Registration QC reports, for subjects whose transforms are computed here."""
    method = config["registration_method"]
    regqc = inputs["spim"].expand(
        rules.registration_qc_report.output.report_html,
        template=config["template"],
        reg=method,
    )
    # same entity order as regqc, so the lists line up per subject
    if cohort_warps:
        # subjects from the cohort build have no registration of their own to check
        built = {
            os.path.basename(f).split("_stain-")[0]
            for f in _glob.glob(
                os.path.join(cohort_warps["warps_dir"], "*_invwarp.nii.gz")
            )
        }
        regqc = [
            qc for qc in regqc if os.path.basename(qc).split("_stain-")[0] not in built
        ]
    if reuse_root and "transforms" in reuse_outputs and not cohort_warps:
        reused = inputs["spim"].expand(
            reuse_path(rules.compose_subject_to_template_warp.output.xfm_composite),
            template=config["template"],
            reg=method,
        )
        legacy = inputs["spim"].expand(
            reuse_path(rules.select_subject_to_template_warp.output.xfm_composite),
            template=config["template"],
        )
        regqc = [
            qc
            for qc, xfm, old in zip(regqc, reused, legacy)
            if not os.path.exists(xfm)
            and not (method == "greedy" and os.path.exists(old))
        ]
    return regqc
