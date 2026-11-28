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

Segmentation and vessel outputs (opt-in: --reuse_outputs segmentation vessels)
depend on settings their paths do not name (e.g. --correction_method for a gmm
mask). Those settings are compared with the prior dataset's
config/snakebids.yml first, and a group of outputs whose settings differ is not
reused at all (see scripts/reuse_check.py). Their QC reports are made as usual:
they are cheap, and show this run's tables.
"""

import glob as _glob
import re as _re
import shlex

from snakemake.io import strip_wildcard_constraints as _strip_constraints

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
# are ~1.5 GB per subject); input paths are absolute, so the links are too.
# OME-Zarr outputs are directories, so they are copied recursively (-L: a prior
# dataset made with --reuse_symlink holds links, copy what they point to).
reuse_cmd = "ln -s" if config.get("reuse_symlink") else "cp"
reuse_tree_cmd = "ln -s" if config.get("reuse_symlink") else "cp -rL"


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
            reg=qc_reg_constraint(),
        localrule: True
        shell:
            f"{reuse_cmd} {{input.xfm_composite}} {{output.xfm_composite}} && "
            f"{reuse_cmd} {{input.xfm_composite_inv}} {{output.xfm_composite_inv}}"

    # (and > the --reg-cohort-template compose rules, in cohort_template.smk)
    ruleorder: reuse_composite_warps > compose_subject_to_template_warp

    if via_template:

        ruleorder: reuse_composite_warps > compose_subject_warp_via_template

    _reuse_method = subject_template_method()
    _reuse_found = []
    _reuse_legacy = 0
    for desc_path, legacy_path in zip(
        ref_spim.expand(
            reuse_path(rules.compose_subject_to_template_warp.output.xfm_composite),
            template=config["template"],
            reg=_reuse_method,
        ),
        ref_spim.expand(
            reuse_path(rules.select_subject_to_template_warp.output.xfm_composite),
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
            "Atlas labels from these warps will inherit that error, so register again "
            "if you need accurate labels."
        )


# --- segmentation and vessel outputs (opt-in) --------------------------------

from reuse_check import (
    REUSE_SETTINGS,
    format_differences,
    load_prior_config,
    setting_differences,
)
from seg_plan import SEG_METHOD_PATTERNS


def _reuse_output(path):
    """Output path of a reuse rule: OME-Zarr stores are directories."""
    path = str(path)
    return directory(path) if path.endswith(".ome.zarr") else path


def _reuse_input(path):
    """The prior dataset's copy of an output path (minus wildcard regexes)."""
    return reuse_path(_strip_constraints(str(path)))


def _alternation(values):
    return "|".join(_re.escape(v) for v in values) or "^$"


def _seg_method_kind(method):
    return next(k for k, p in SEG_METHOD_PATTERNS.items() if p.match(method))


def _reuse_settings_ok(group, what):
    """Whether a group's settings match the prior dataset's (says why if not)."""
    if _prior_config is None:
        return False
    diffs = setting_differences(_prior_config, config, REUSE_SETTINGS[group])
    if diffs:
        print(
            f"Not reusing {what} from {reuse_root}, settings differ: "
            f"{format_differences(diffs)}. They will be computed."
        )
    return not diffs


# files written together with a seg method's mask, by kind of method: a mask is
# only reused when these are in the prior dataset too
_seg_mask_siblings = {
    "otsu": ["thresholds.png"],
    "gmm": ["thresholds.png"],
    "dog": ["thresholds.png", "poolmask.nii.gz"],
}


def _with_siblings(mask_path, siblings):
    base = str(mask_path).rsplit("_mask.", 1)[0]
    return [mask_path] + [f"{base}_{s}" for s in siblings]


def _print_found(what, path_sets, unit):
    """How many of ``path_sets`` (files reused together) the prior dataset has."""
    n_found = sum(
        all(os.path.exists(reuse_path(p)) for p in paths) for paths in path_sets
    )
    print(
        f"Reusing {what} from {reuse_root} for {n_found} of {len(path_sets)} "
        f"{unit}; the rest will be computed"
    )


# reuse rules to define: (name, producing rule's name, wildcard constraints)
_reuse_specs = []
_prior_config = None

if reuse_root and {"segmentation", "vessels"} & set(reuse_outputs):
    _prior_config = load_prior_config(reuse_root)
    if _prior_config is None:
        print(
            f"Not reusing segmentation or vessel outputs from {reuse_root}: it has "
            "no config/snakebids.yml to check their settings against"
        )


def _seg_mask_paths(where):
    return seg_expand(
        bids_oz_in(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
        kind="mask",
        where=where,
    )


if _prior_config is not None and "segmentation" in reuse_outputs:
    _kinds = {m: _seg_method_kind(m) for m in seg_plan.all_seg_methods()}
    # methods whose masks (and SDTs) are reused
    _reused_methods = []

    def _stains_of_kind(kind):
        return [
            s for m, k in _kinds.items() if k == kind for s in seg_plan.stains_using(m)
        ]

    _threshold_kinds = [k for k in ("th", "otsu", "gmm") if k in _kinds.values()]
    if _threshold_kinds and _reuse_settings_ok(
        "threshold", "th/otsu/gmm segmentation masks"
    ):
        for kind, producer in [
            ("th", "threshold"),
            ("otsu", "multiotsu"),
            ("gmm", "gmmthresh"),
        ]:
            if kind in _threshold_kinds:
                _reuse_specs.append(
                    (
                        f"reuse_{producer}",
                        producer,
                        {"stain": _alternation(_stains_of_kind(kind))},
                    )
                )
        _reused_methods += [m for m, k in _kinds.items() if k in _threshold_kinds]

    if "dog" in _kinds.values() and _reuse_settings_ok("dog", "dog segmentation masks"):
        _reuse_specs.append(
            (
                "reuse_dogthresh",
                "dogthresh",
                {"stain": _alternation(_stains_of_kind("dog"))},
            )
        )
        _reused_methods += [m for m, k in _kinds.items() if k == "dog"]

    if lantern_stains and _reuse_settings_ok(
        "lantern_probseg", "LANTERN vote maps (probseg)"
    ):
        _reuse_specs.append(
            (
                "reuse_lantern_probseg",
                "run_lantern_plaques",
                {"stain": _alternation(lantern_stains)},
            )
        )
        _print_found(
            "LANTERN vote maps",
            [
                [p]
                for p in inputs["spim"].expand(
                    rules.binarize_lantern_plaques.input.probseg, stain=lantern_stains
                )
            ],
            "scan(s)",
        )
        if _reuse_settings_ok("lantern_mask", "LANTERN plaque masks"):
            _reuse_specs.append(
                (
                    "reuse_lantern_mask",
                    "binarize_lantern_plaques",
                    {"stain": _alternation(lantern_stains)},
                )
            )
            _reused_methods += [m for m, k in _kinds.items() if k == "lantern"]

    if _reused_methods:
        _print_found(
            f"segmentation masks ({', '.join(_reused_methods)})",
            [
                _with_siblings(p, _seg_mask_siblings.get(_kinds[m], []))
                for m in _reused_methods
                for p in _seg_mask_paths(lambda method: method == m)
            ],
            "mask(s)",
        )
        # signed distance transforms of these masks (sdt_stains)
        _reuse_specs.append(
            (
                "reuse_seg_dist",
                "signed_distance_transform",
                {
                    "datatype": "seg",
                    "stain": _alternation(list(seg_plan.stains)),
                    "desc": _alternation(_reused_methods),
                },
            )
        )

if (
    _prior_config is not None
    and "vessels" in reuse_outputs
    and do_vessels
    and _reuse_settings_ok("vessels", "vessel outputs")
):
    _vessel_constraints = {
        "stain": _re.escape(stain_for_vessels),
        "desc": _re.escape(config["vessel_seg_method"]),
    }
    for producer in [
        "run_vesselfm",
        "skeletonize_vessels_mask",
        "vessel_skeleton_graph",
        "vessel_graph_to_nodes_edges",
    ]:
        _reuse_specs.append((f"reuse_{producer}", producer, _vessel_constraints))
    _reuse_specs.append(
        (
            "reuse_vessel_dist",
            "signed_distance_transform",
            {"datatype": "vessels", **_vessel_constraints},
        )
    )
    _print_found(
        "vessel masks",
        [
            [p]
            for p in inputs["spim"].expand(
                str(rules.run_vesselfm.output.mask),
                stain=stain_for_vessels,
                level=config["segmentation_level"],
            )
        ],
        "scan(s)",
    )


# One rule per producing rule, copying (or linking) all of its outputs from the
# prior dataset, so the producing rule is never needed for one output while the
# reuse rule makes another. Where the outputs are not all in the prior dataset,
# the input is missing and Snakemake falls back to the producing rule.
for _name, _producer_name, _constraints in _reuse_specs:
    _producer = getattr(rules, _producer_name)

    rule:
        name:
            _name
        input:
            **{k: _reuse_input(v) for k, v in _producer.output.items()},
        output:
            # the SDT rule's inline {datatype,seg|vessels} would override the
            # constraint here; elsewhere keep the producer's inline regexes
            **{
                k: _reuse_output(
                    _strip_constraints(str(v)) if "datatype" in _constraints else v
                )
                for k, v in _producer.output.items()
            },
        wildcard_constraints:
            **_constraints,
        localrule: True
        params:
            # a function with the keys bound now: rules defined in a loop
            # evaluate their shell string when the job runs, after the loop
            cmd=lambda wildcards, input, output, keys=list(
                _producer.output.keys()
            ): " && ".join(
                f"{reuse_tree_cmd} {shlex.quote(input[k])} {shlex.quote(output[k])}"
                for k in keys
            ),
        shell:
            "{params.cmd}"

    workflow.ruleorder(_name, _producer_name)


def get_regqc_targets(wildcards):
    """Registration QC reports, for subjects whose transforms are computed here.

    Only reference scans are registered to the template; additional scans
    (--aux_spim) get the reference scan resampled into them instead.
    With --via-template, the reports are of the registrations to the via template.
    """
    method = (
        config["registration_method"] if via_template else subject_template_method()
    )
    regqc = ref_spim.expand(
        rules.registration_qc_report.output.report_html,
        template=reg_template,
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
        reused = ref_spim.expand(
            reuse_path(rules.compose_subject_to_template_warp.output.xfm_composite),
            template=reg_template,
            reg=method,
        )
        legacy = ref_spim.expand(
            reuse_path(rules.select_subject_to_template_warp.output.xfm_composite),
            template=reg_template,
        )
        regqc = [
            qc
            for qc, xfm, old in zip(regqc, reused, legacy)
            if not os.path.exists(xfm)
            and not (method == "greedy" and os.path.exists(old))
        ]
    if reg_via_cohort and regqc and not prior_cohort:
        # the cohort template's registration to --template, built in this run
        regqc.append(rules.cohort_template_regqc.output.report_html)
    return regqc
