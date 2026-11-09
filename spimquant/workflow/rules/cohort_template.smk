"""
Cohort template building workflow for SPIMquant.

Builds a study-specific (cohort) SPIM template from subjects affine-registered to
--template, then registers the original template to the new one so its atlases can
be carried over into cohort template space.

Key workflow stages:
1. Resample each subject's per-stain images to --template space with an affine
   (plus images from any --cohort-template-extra-inputs SPIMquant output directories)
2. Robust intensity normalization of the subject images
3. Multivariate template building with FireANTs (all stains share one space)
   and a desc-rawintensity template per stain: the subjects' un-normalized images
   warped with the final-iteration warps and averaged
4. Registration of the original template (moving) to the cohort template
   (fixed), since the cohort template drifts away from the original template
   space during building, with --cohort-template-reg-method. The
   cohort_template_tune_reg analysis level compares methods side by side.
5. Warping of the original atlases, brain mask and anat into cohort template space
6. Export: template_description.json, a templates: config snippet for SPIMquant,
   and optional upload to a Zenodo deposition

Invoke with:
``spimquant <bids_dir> <out_dir> cohort_template`` and then
``spimquant <bids_dir> <out_dir> cohort_template_export``
"""

import glob as _glob

cohort_name = config["cohort_template_name"]
orig_template = config["template"]
cohort_level = (
    config["cohort_template_level"]
    if config["cohort_template_level"] is not None
    else config["registration_level"]
)
cohort_iterations = config["cohort_template_iterations"]
cohort_n_gpus = config["cohort_template_n_gpus"]
cohort_zenodo_deposition = config["cohort_template_zenodo_deposition"]
cohort_reg_method = config["cohort_template_reg_method"]
cohort_reg_tune_methods = config["cohort_template_reg_tune_methods"]

# final outputs in cohort_root/tpl-{name}/, intermediate files in cohort_work
cohort_root = os.path.join(root, "cohort_template")
cohort_work = os.path.join(work, "cohort_template", f"work-{cohort_name}")

# the registration stain is channel 0 of the multivariate template: it is the
# fixed image for the original-to-cohort template registration
_cohort_stains_req = config["cohort_template_stains"] or stains
cohort_stains = [stain_for_reg] if stain_for_reg in _cohort_stains_req else []
cohort_stains += [s for s in _cohort_stains_req if s not in cohort_stains]
cohort_reg_stain = cohort_stains[0]

# atlases carried over from the original template (all of them, so the cohort
# template supports every --atlas_segs choice the original did)
cohort_atlas_segs = list(template.get("atlases", {}).keys())
cohort_has_mask = "mask" in template

# snakebids paths always use the bids() entity order, so the filename tail is
# the same in other SPIMquant output directories
_cohort_input_pattern = bids(
    root=root,
    datatype="micr",
    stain="{stain}",
    level="{level}",
    space="{template}",
    desc="affine",
    suffix="SPIM.nii.gz",
    **inputs["spim"].wildcards,
)


def _cohort_row_name(path, stain):
    """Subject identifier: the filename entities that precede stain-."""
    return os.path.basename(path).split(f"_stain-{stain}")[0]


def _get_cohort_rows():
    """One row per subject: name and per-stain image paths (cohort_stains order)."""
    rows = []
    per_stain = [
        inputs["spim"].expand(
            _cohort_input_pattern,
            stain=s,
            level=cohort_level,
            template=orig_template,
        )
        for s in cohort_stains
    ]
    for images in zip(*per_stain):
        rows.append(
            {
                "name": _cohort_row_name(images[0], cohort_stains[0]),
                "images": list(images),
            }
        )

    for extra_dir in config["cohort_template_extra_inputs"] or []:
        by_name = {}
        for s in cohort_stains:
            tail = os.path.basename(
                bids(
                    stain=s,
                    level=cohort_level,
                    space=orig_template,
                    suffix="SPIM.nii.gz",
                )
            )
            for path in sorted(
                _glob.glob(
                    os.path.join(extra_dir, "sub-*", "**", f"*_{tail}"), recursive=True
                )
            ):
                by_name.setdefault(_cohort_row_name(path, s), {})[s] = path
        if not by_name:
            raise ValueError(
                f"--cohort-template-extra-inputs: no space-{orig_template} "
                f"level-{cohort_level} SPIM images found in {extra_dir}"
            )
        for name, images in sorted(by_name.items()):
            missing = [s for s in cohort_stains if s not in images]
            if missing:
                print(
                    f"[cohort_template] skipping {name} from {extra_dir}: "
                    f"missing stain(s) {missing}"
                )
                continue
            rows.append({"name": name, "images": [images[s] for s in cohort_stains]})

    # subject IDs can repeat across datasets
    seen = {}
    for row in rows:
        n = seen.get(row["name"], 0)
        seen[row["name"]] = n + 1
        if n > 0:
            row["name"] = f"{row['name']}-dup{n}"
    return rows


if config.get("analysis_level") in (
    "cohort_template",
    "cohort_template_tune_reg",
    "cohort_template_export",
):
    for _method in [cohort_reg_method, *cohort_reg_tune_methods]:
        if _method not in reg_presets:
            raise ValueError(
                f"Unknown registration method {repr(_method)}, choose from "
                f"{sorted(reg_presets)}"
            )
    if not cohort_name.isalnum():
        raise ValueError(
            f"--cohort-template-name must be alphanumeric, got {repr(cohort_name)}"
        )
    if cohort_name in config["templates"]:
        raise ValueError(
            f"--cohort-template-name {repr(cohort_name)} clashes with an existing template"
        )
    _unknown = set(_cohort_stains_req) - set(stains)
    if _unknown:
        raise ValueError(
            f"--cohort-template-stains {sorted(_unknown)} not found in stains {stains}"
        )
    cohort_rows = _get_cohort_rows()
    if len(cohort_rows) < 2:
        raise ValueError(
            f"Need at least 2 subjects to build a cohort template, found {len(cohort_rows)}"
        )
    print(
        f"[cohort_template] building tpl-{cohort_name} from {len(cohort_rows)} subjects, "
        f"stains {cohort_stains}, in {orig_template} space at level {cohort_level}"
    )
else:
    cohort_rows = []


def bids_cohort(**kwargs):
    """bids() path for a final cohort template file (cohort_root/tpl-{name}/)."""
    return bids(root=cohort_root, template=cohort_name, **kwargs)


def bids_orig(root=root, **kwargs):
    """bids() path for an original template file (imported ones by default)."""
    return bids(root=root, template=orig_template, **kwargs)


cohort_description = os.path.join(
    cohort_root, f"tpl-{cohort_name}", "template_description.json"
)
cohort_config_snippet = os.path.join(cohort_root, f"tpl-{cohort_name}_spimquant.yml")
cohort_buildinfo = os.path.join(cohort_root, f"tpl-{cohort_name}_buildinfo.tar.gz")


rule all_cohort_template:
    input:
        expand(bids_cohort(stain="{stain}", suffix="SPIM.nii.gz"), stain=cohort_stains),
        expand(
            bids_cohort(desc="rawintensity", stain="{stain}", suffix="SPIM.nii.gz"),
            stain=cohort_stains,
        ),
        bids_cohort(suffix="anat.nii.gz"),
        bids_cohort(desc="brain", suffix="mask.nii.gz") if cohort_has_mask else [],
        expand(
            bids_cohort(seg="{seg}", suffix="dseg.{ext}"),
            seg=cohort_atlas_segs,
            ext=["nii.gz", "tsv"],
        ),
        bids_cohort(
            from_=orig_template,
            to=cohort_name,
            desc=cohort_reg_method,
            suffix="xfm.nii.gz",
        ),
        bids_cohort(
            from_=cohort_name,
            to=orig_template,
            desc=cohort_reg_method,
            suffix="xfm.nii.gz",
        ),
        bids_cohort(
            stain=cohort_reg_stain,
            from_=orig_template,
            to=cohort_name,
            desc=cohort_reg_method,
            suffix="regqc.html",
        ),
        # cohort template resampled back into the original template space
        expand(
            bids_orig(
                root=cohort_root,
                desc=cohort_name,
                stain="{stain}",
                suffix="SPIM.nii.gz",
            ),
            stain=cohort_stains,
        ),


rule cohort_template_affine_to_itk:
    """Convert a subject's initial affine to the original template to ITK format."""
    input:
        xfm_ras=rules.init_affine_reg.output.xfm_ras,
    output:
        xfm_itk=temp(_reg_xfm("xfm.txt", desc="initaffine", type_="itk")),
    group:
        "cohort_affine"
    conda:
        "../envs/c3d.yaml"
    threads: 1
    resources:
        mem_mb=1500,
        runtime=10,
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


rule cohort_template_affine_warp:
    """Composite warps for a subject's initial affine to the original template.

    Building needs only an affine, since the template gets its own shape, and the
    original template is registered to it afterwards. Participant runs with the
    cohort template chain these with the warps saved by the build.
    """
    input:
        ref=rules.import_template_anat.output.anat,
        subject=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=config["registration_level"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        xfm_itk=rules.cohort_template_affine_to_itk.output.xfm_itk,
    output:
        xfm_composite=_reg_xfm("xfm.nii.gz", desc="affine"),
        xfm_composite_inv=_reg_xfm(
            "xfm.nii.gz", from_="{template}", to="subject", desc="affine"
        ),
    group:
        "cohort_affine"
    conda:
        "../envs/ants.yaml"
    threads: 4
    resources:
        mem_mb=10000,
        runtime=15,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite},1] "
        " -r {input.ref} -t {input.xfm_itk} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite_inv},1] "
        " -r {input.subject} -t [{input.xfm_itk},1]"


use rule deform_spim_nii_to_template_nii as cohort_template_affine_input with:
    input:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        ref=rules.import_template_anat.output.anat,
        xfm_composite=rules.cohort_template_affine_warp.output.xfm_composite,
    output:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            space="{template}",
            desc="affine",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    group:
        "cohort_affine"


rule cohort_template_inputs:
    """Normalize subject images and write the FireANTs image list and config.

    Images are optionally resampled to --cohort-template-resample-mm. Each image
    is clipped to robust percentiles of its nonzero voxels and rescaled to [0,1]. The image list has one comma-separated row per subject
    (one entry per stain), which FireANTs builds as one multivariate template.
    """
    input:
        images=lambda wildcards: [img for row in cohort_rows for img in row["images"]],
    output:
        images=temp(directory(os.path.join(cohort_work, "inputs"))),
        # temp: only build_cohort_template reads these, and it keeps copies in
        # {build_dir}/buildinfo/ for provenance, so losing them never forces a
        # rebuild of the template
        image_list=temp(os.path.join(cohort_work, "images.txt")),
        fireants_config=temp(os.path.join(cohort_work, "config", "cohort.yaml")),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=lambda wildcards, input: 10 + len(input.images),
    params:
        rows=cohort_rows,
        stains=cohort_stains,
        clip_percentiles=config["cohort_template"]["clip_percentiles"],
        resample_mm=config["cohort_template_resample_mm"],
        fireants=lambda wildcards, output: {
            **config["cohort_template"]["fireants"],
            "template_iterations": cohort_iterations,
            "image_list_file": os.path.abspath(output.image_list),
            "save_dir": os.path.abspath(
                os.path.join(cohort_work, "fireants", "buildtemplate")
            ),
        },
    script:
        "../scripts/cohort_template_inputs.py"


rule build_cohort_template:
    """Build the multivariate cohort template with FireANTs (GPU).

    Also writes each subject's final-iteration warp to {build_dir}/warps/, and
    copies the image list and FireANTs config to {build_dir}/buildinfo/ so that
    downstream (export) rules never depend on the temp inputs.
    """
    input:
        images=rules.cohort_template_inputs.output.images,
        image_list=rules.cohort_template_inputs.output.image_list,
        fireants_config=rules.cohort_template_inputs.output.fireants_config,
    output:
        build_dir=directory(os.path.join(cohort_work, "fireants", "buildtemplate")),
    log:
        os.path.join(cohort_work, "logs", "build_cohort_template.log"),
    threads: 8 * cohort_n_gpus
    resources:
        gpu=cohort_n_gpus,
        # one task that launches one process per GPU via torchrun (see plaques.smk)
        tasks_per_gpu=0,
        cpus_per_gpu=8,
        mem_mb=64000,
        # rough estimate: ~30s per subject per iteration, split across GPUs
        runtime=int(15 + 0.5 * len(cohort_rows) * cohort_iterations / cohort_n_gpus),
    params:
        config_dir=lambda wildcards, input: os.path.abspath(
            os.path.dirname(input.fireants_config)
        ),
        config_name="cohort",
        hydra_dir=lambda wildcards, output: os.path.abspath(
            os.path.join(output.build_dir, "..", "hydra")
        ),
    shell:
        # FireANTs is only in the gpu pixi environment
        "python -c 'import fireants' 2>/dev/null || {{ echo 'FireANTs is not "
        "installed: run SPIMquant from the gpu environment, e.g. pixi run -e gpu "
        "spimquant ...' >&2; exit 1; }} && "
        # full-res deformable steps can run close to the GPU memory limit, and
        # fragmentation otherwise causes OOM after a few subjects
        "PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True "
        "torchrun --standalone --nproc_per_node {resources.gpu} "
        # script path in the shell (not a param), so moving or reinstalling
        # SPIMquant does not count as a params change and rebuild the template
        " {workflow.basedir}/scripts/cohort_template_build.py "
        " --config-dir {params.config_dir} --config-name {params.config_name} "
        " hydra.run.dir={params.hydra_dir} &> {log} && "
        "mkdir -p {output.build_dir}/buildinfo && "
        "cp {input.image_list} {input.fireants_config} {output.build_dir}/buildinfo/"


rule cohort_template_stain:
    """Extract the final-iteration template for one stain (channel)."""
    input:
        build_dir=rules.build_cohort_template.output.build_dir,
    output:
        spim=bids_cohort(stain="{stain}", suffix="SPIM.nii.gz"),
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    params:
        template_nii=lambda wildcards, input: os.path.join(
            input.build_dir,
            f"template_{cohort_iterations-1}"
            + (
                f"_ch-{cohort_stains.index(wildcards.stain)}"
                if len(cohort_stains) > 1
                else ""
            )
            + ".nii.gz",
        ),
    shell:
        "if [ ! -f {params.template_nii} ]; then "
        " echo 'Missing {params.template_nii}: multi-stain cohort templates need a "
        "FireANTs build_template that saves every channel (template_<epoch>_ch-<c>)' >&2;"
        " exit 1; fi; "
        "cp {params.template_nii} {output.spim}"


_cohort_rows_by_name = {row["name"]: row for row in cohort_rows}
_cohort_raw_warped = os.path.join(
    cohort_work, "rawintensity", "{row}_stain-{stain}.nii.gz"
)


rule cohort_template_warp_rawintensity:
    """Warp a subject's un-normalized image with its final-iteration warp.

    The image is the input to cohort_template_inputs, before intensity clipping
    and rescaling, so intensities stay comparable across subjects.
    """
    input:
        build_dir=rules.build_cohort_template.output.build_dir,
        img=lambda wildcards: _cohort_rows_by_name[wildcards.row]["images"][
            cohort_stains.index(wildcards.stain)
        ],
        ref=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
    output:
        temp(_cohort_raw_warped),
    wildcard_constraints:
        row="[^/]+",
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        # inside the build directory, so it cannot be declared as an input
        warp=lambda wildcards, input: os.path.join(
            input.build_dir,
            "warps",
            f"{wildcards.row}_stain-{cohort_stains[0]}_warp.nii.gz",
        ),
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -v -n Linear -u float "
        " -i {input.img} -o {output} -r {input.ref} -t {params.warp}"


rule cohort_template_rawintensity:
    """Average the warped un-normalized images into a desc-rawintensity template.

    A plain voxelwise mean, without the Laplacian sharpening and min-max
    normalization FireANTs applies to its template.
    """
    input:
        images=lambda wildcards: expand(
            _cohort_raw_warped,
            row=[row["name"] for row in cohort_rows],
            stain=wildcards.stain,
        ),
    output:
        bids_cohort(desc="rawintensity", stain="{stain}", suffix="SPIM.nii.gz"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=lambda wildcards, input: 10 + len(input.images),
    script:
        "../scripts/cohort_template_rawintensity.py"


_cohort_moving = bids_orig(suffix="anat.nii.gz")
_cohort_fixed = bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz")
_cohort_affine_ras = os.path.join(cohort_work, "xfm", "{reg}_affine_ras_xfm.txt")
_cohort_affine_itk = os.path.join(cohort_work, "xfm", "{reg}_affine_itk_xfm.txt")
_cohort_warp = os.path.join(cohort_work, "xfm", "{reg}_warp.nii.gz")
_cohort_invwarp = os.path.join(cohort_work, "xfm", "{reg}_invwarp.nii.gz")
_cohort_warped_affine = os.path.join(cohort_work, "xfm", "{reg}_affinewarped.nii.gz")
_cohort_warped = os.path.join(cohort_work, "xfm", "{reg}_deformwarped.nii.gz")


def _cohort_extra_fixed(wildcards):
    """Other stains of the cohort template, for multichannel registration methods."""
    if not reg_presets[wildcards.reg].get("multichannel"):
        return []
    return [
        bids_cohort(stain=stain, suffix="SPIM.nii.gz")
        for stain in cohort_stains
        if stain != cohort_reg_stain
    ]


def _cohort_extra_pairs(wildcards):
    """greedy image pairs for the other stains, each against the same moving image."""
    return " ".join(
        f"-i {fixed} {_cohort_moving}" for fixed in _cohort_extra_fixed(wildcards)
    )


def _cohort_greedy_mem_mb(wildcards):
    return 32000 if reg_presets[wildcards.reg].get("multichannel") else 16000


_cohort_reg_log = os.path.join(cohort_work, "logs", "{reg}_{step}.log")


# The images start out roughly aligned (subjects were affine-registered to the
# original template), so registration is initialized with identity
use rule affine_reg as cohort_template_affine_reg with:
    input:
        fixed=_cohort_fixed,
        moving=_cohort_moving,
        extra_fixed=_cohort_extra_fixed,
    output:
        xfm_ras=temp(_cohort_affine_ras),
        warped=temp(_cohort_warped_affine),
    log:
        _cohort_reg_log.replace("{step}", "affine_reg"),
    group:
        "cohort_reg"
    resources:
        mem_mb=_cohort_greedy_mem_mb,
        runtime=15,
    params:
        affine=lambda wildcards: reg_presets[wildcards.reg]["affine"],
        init="-ia-identity",
        metric=lambda wildcards: reg_presets[wildcards.reg]["metric"],
        iters=lambda wildcards: reg_iters(reg_presets[wildcards.reg], "affine_iters"),
        identity=reg_identity_ras,
        extra_pairs=_cohort_extra_pairs,


use rule convert_ras_to_itk as cohort_template_ras_to_itk with:
    input:
        xfm_ras=_cohort_affine_ras,
    output:
        xfm_itk=temp(_cohort_affine_itk),
    group:
        "cohort_reg"


use rule deform_reg as cohort_template_deform_reg with:
    input:
        fixed=_cohort_fixed,
        moving=_cohort_moving,
        extra_fixed=_cohort_extra_fixed,
        xfm_ras=_cohort_affine_ras,
    output:
        warp=temp(_cohort_warp),
        invwarp=temp(_cohort_invwarp),
        warped=temp(_cohort_warped),
    log:
        _cohort_reg_log.replace("{step}", "deform_reg"),
    group:
        "cohort_reg"
    resources:
        mem_mb=_cohort_greedy_mem_mb,
        runtime=5 if config["sloppy"] else 60,
    params:
        iters=lambda wildcards: reg_iters(reg_presets[wildcards.reg], "deform_iters"),
        metric=lambda wildcards: reg_presets[wildcards.reg]["metric"],
        sigma1=lambda wildcards: reg_presets[wildcards.reg]["sigma1"],
        sigma2=lambda wildcards: reg_presets[wildcards.reg]["sigma2"],
        extra_pairs=_cohort_extra_pairs,


use rule reg_ants as cohort_template_reg_ants with:
    input:
        fixed=_cohort_fixed,
        moving=_cohort_moving,
        extra_fixed=_cohort_extra_fixed,
    output:
        xfm_itk=temp(_cohort_affine_itk),
        warp=temp(_cohort_warp),
        invwarp=temp(_cohort_invwarp),
        warped_affine=temp(_cohort_warped_affine),
        warped=temp(_cohort_warped),
    log:
        _cohort_reg_log.replace("{step}", "reg_ants"),
    params:
        args=lambda wildcards, input: ants_reg_args(
            reg_presets[wildcards.reg],
            input.fixed,
            input.moving,
            center_init=False,
            extra_fixed=_cohort_extra_fixed(wildcards),
        ),
        identity=reg_identity_itk,


use rule reg_fireants as cohort_template_reg_fireants with:
    input:
        fixed=_cohort_fixed,
        moving=_cohort_moving,
        extra_fixed=_cohort_extra_fixed,
    output:
        xfm_itk=temp(_cohort_affine_itk),
        warp=temp(_cohort_warp),
        invwarp=temp(_cohort_invwarp),
        warped_affine=temp(_cohort_warped_affine),
        warped=temp(_cohort_warped),
    params:
        preset=lambda wildcards: {
            **reg_presets[wildcards.reg],
            "affine_iters": reg_iters(reg_presets[wildcards.reg], "affine_iters"),
            "deform_iters": reg_iters(reg_presets[wildcards.reg], "deform_iters"),
        },
        init_centers=False,
        extra_fixed=_cohort_extra_fixed,
        extra_moving=[],


rule cohort_template_compose_warp:
    """Composite warps between the original and cohort template spaces.

    from-{orig}_to-{cohort} resamples original template images into cohort
    template space (and maps cohort-space points to the original space);
    from-{cohort}_to-{orig} is its inverse.
    """
    input:
        cohort_ref=_cohort_fixed,
        orig_ref=bids_orig(suffix="anat.nii.gz"),
        xfm_itk=_cohort_affine_itk,
        warp=_cohort_warp,
        invwarp=_cohort_invwarp,
    output:
        xfm_composite=bids_cohort(
            from_=orig_template, to=cohort_name, desc="{reg}", suffix="xfm.nii.gz"
        ),
        xfm_composite_inv=bids_cohort(
            from_=cohort_name, to=orig_template, desc="{reg}", suffix="xfm.nii.gz"
        ),
    wildcard_constraints:
        reg=reg_constraint(),
    conda:
        "../envs/ants.yaml"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=15,
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite},1] "
        " -r {input.cohort_ref} -t {input.warp} {input.xfm_itk} && "
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "antsApplyTransforms -d 3 -o [{output.xfm_composite_inv},1] "
        " -r {input.orig_ref} -t {input.invwarp} [{input.xfm_itk},1]"


# the warps of the method selected with --cohort-template-reg-method
_cohort_xfm = bids_cohort(
    from_=orig_template, to=cohort_name, desc=cohort_reg_method, suffix="xfm.nii.gz"
)
_cohort_xfm_inv = bids_cohort(
    from_=cohort_name, to=orig_template, desc=cohort_reg_method, suffix="xfm.nii.gz"
)


_cohort_warp_shell = (
    "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
    "antsApplyTransforms -d 3 -v -n {params.interp} "
    " -i {input.img} -o {output} -r {input.ref} -t {input.xfm}"
)


rule cohort_template_warp_anat:
    input:
        img=bids_orig(suffix="anat.nii.gz"),
        ref=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        xfm=_cohort_xfm,
    output:
        bids_cohort(suffix="anat.nii.gz"),
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        interp="Linear",
    shell:
        _cohort_warp_shell


rule cohort_template_warp_mask:
    input:
        img=bids_orig(desc="brain", suffix="mask.nii.gz"),
        ref=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        xfm=_cohort_xfm,
    output:
        bids_cohort(desc="brain", suffix="mask.nii.gz"),
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        interp="NearestNeighbor",
    shell:
        _cohort_warp_shell


rule cohort_template_warp_dseg:
    """Warp an original template atlas into cohort template space."""
    input:
        img=bids_orig(seg="{seg}", suffix="dseg.nii.gz"),
        ref=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        xfm=_cohort_xfm,
    output:
        bids_cohort(seg="{seg}", suffix="dseg.nii.gz"),
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        interp="NearestNeighbor",
    shell:
        _cohort_warp_shell


rule cohort_template_copy_dseg_tsv:
    input:
        bids_orig(seg="{seg}", suffix="dseg.tsv"),
    output:
        bids_cohort(seg="{seg}", suffix="dseg.tsv"),
    localrule: True
    shell:
        "cp {input} {output}"


rule cohort_template_to_orig_space:
    """Resample a cohort stain template into the original template space.

    Gives stain templates that are aligned to the original template's atlases.
    """
    input:
        img=bids_cohort(stain="{stain}", suffix="SPIM.nii.gz"),
        ref=bids_orig(suffix="anat.nii.gz"),
        xfm=_cohort_xfm_inv,
    output:
        bids_orig(
            root=cohort_root, desc=cohort_name, stain="{stain}", suffix="SPIM.nii.gz"
        ),
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        interp="Linear",
    shell:
        _cohort_warp_shell


_cohort_report_stains = [cohort_reg_stain] + [
    stain for stain in cohort_stains if stain != cohort_reg_stain
]
_cohort_report_extra_fixed = [
    bids_cohort(stain=stain, suffix="SPIM.nii.gz")
    for stain in _cohort_report_stains[1:]
]


rule cohort_template_regqc:
    """QC report for the original-to-cohort template registration."""
    input:
        template=_cohort_fixed,
        extra_fixed=_cohort_report_extra_fixed,
        warped_affine=_cohort_warped_affine.format(reg=cohort_reg_method),
        warped_deform=_cohort_warped.format(reg=cohort_reg_method),
        warp=_cohort_xfm,
        dseg=bids_cohort(seg=list(atlas_segs)[0], suffix="dseg.nii.gz"),
    output:
        report_html=bids_cohort(
            stain=cohort_reg_stain,
            from_=orig_template,
            to=cohort_name,
            desc=cohort_reg_method,
            suffix="regqc.html",
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        subject=f"tpl-{orig_template}",
        template=f"tpl-{cohort_name}",
        stains=_cohort_report_stains,
        method=cohort_reg_method,
        parameters=reg_note(reg_presets[cohort_reg_method]),
        labels_warped=True,
    script:
        "../scripts/reg_qc_report.py"


rule cohort_template_tune_warp_dseg:
    """Warp an original template atlas with a method compared by the tune report."""
    input:
        img=bids_orig(seg="{seg}", suffix="dseg.nii.gz"),
        ref=_cohort_fixed,
        xfm=bids_cohort(
            from_=orig_template, to=cohort_name, desc="{reg}", suffix="xfm.nii.gz"
        ),
    output:
        temp(os.path.join(cohort_work, "tune", "{reg}_seg-{seg}_dseg.nii.gz")),
    wildcard_constraints:
        reg=reg_constraint(),
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        interp="NearestNeighbor",
    shell:
        _cohort_warp_shell


rule cohort_template_tune_warp_anat:
    """Warp the original template with a method compared by the tune report.

    From the composite warp, which is kept, so that the report can be written again
    without repeating the registrations.
    """
    input:
        img=bids_orig(suffix="anat.nii.gz"),
        ref=_cohort_fixed,
        xfm=bids_cohort(
            from_=orig_template, to=cohort_name, desc="{reg}", suffix="xfm.nii.gz"
        ),
    output:
        temp(os.path.join(cohort_work, "tune", "{reg}_anat.nii.gz")),
    wildcard_constraints:
        reg=reg_constraint(),
    conda:
        "../envs/ants.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        interp="Linear",
    shell:
        _cohort_warp_shell


rule cohort_template_tune_report:
    """Compare methods for registering the original template to the cohort template."""
    input:
        fixed=_cohort_fixed,
        extra_fixed=_cohort_report_extra_fixed,
        warped=expand(
            rules.cohort_template_tune_warp_anat.output[0], reg=cohort_reg_tune_methods
        ),
        warps=expand(
            bids_cohort(
                from_=orig_template, to=cohort_name, desc="{reg}", suffix="xfm.nii.gz"
            ),
            reg=cohort_reg_tune_methods,
        ),
        labels=expand(
            rules.cohort_template_tune_warp_dseg.output[0],
            reg=cohort_reg_tune_methods,
            seg=list(atlas_segs)[0],
        ),
    output:
        report_html=bids_cohort(
            stain=cohort_reg_stain,
            from_=orig_template,
            to=cohort_name,
            desc="tune",
            suffix="regqc.html",
        ),
    threads: 1
    resources:
        mem_mb=32000,
        runtime=60,
    params:
        subject=f"tpl-{orig_template}",
        template=f"tpl-{cohort_name}",
        stains=_cohort_report_stains,
        methods=cohort_reg_tune_methods,
        selected=cohort_reg_method,
        default_layer="fixed",
        notes={m: reg_note(reg_presets[m]) for m in reg_presets},
    script:
        "../scripts/reg_tune_report.py"


rule all_cohort_template_reg_tune:
    input:
        rules.cohort_template_tune_report.output.report_html,


# --- export: prepare the cohort template for inclusion in SPIMquant ---

# every file a SPIMquant templates: entry for the cohort template refers to
cohort_export_files = (
    [bids_cohort(suffix="anat.nii.gz")]
    + ([bids_cohort(desc="brain", suffix="mask.nii.gz")] if cohort_has_mask else [])
    + [bids_cohort(stain=s, suffix="SPIM.nii.gz") for s in cohort_stains]
    + [
        bids_cohort(desc="rawintensity", stain=s, suffix="SPIM.nii.gz")
        for s in cohort_stains
    ]
    + [
        bids_cohort(seg=seg, suffix=f"dseg.{ext}")
        for seg in cohort_atlas_segs
        for ext in ["nii.gz", "tsv"]
    ]
    + [_cohort_xfm, _cohort_xfm_inv]
)


cohort_upload_files = cohort_export_files + [cohort_description, cohort_buildinfo]


rule all_cohort_template_export:
    input:
        rules.all_cohort_template.input,
        cohort_description,
        cohort_config_snippet,
        cohort_buildinfo,
        (
            [
                os.path.join(cohort_work, "zenodo", f"{os.path.basename(f)}.uploaded")
                for f in cohort_upload_files
            ]
            if cohort_zenodo_deposition is not None
            else []
        ),


rule cohort_template_buildinfo:
    """Archive of the inputs used to build the template, for provenance."""
    input:
        # logs cannot be inputs; these outputs guarantee the logs were written
        # (build_dir also holds buildinfo/, the image list and FireANTs config)
        build_dir=rules.build_cohort_template.output.build_dir,
        warp=_cohort_xfm,
    output:
        cohort_buildinfo,
    localrule: True
    params:
        buildinfo=lambda wildcards, input: [
            os.path.join(input.build_dir, "buildinfo", f)
            for f in ("images.txt", "cohort.yaml")
        ],
        logs=[rules.build_cohort_template.log[0]]
        + {
            "greedy": [
                rules.cohort_template_affine_reg.log[0].format(reg=cohort_reg_method),
                rules.cohort_template_deform_reg.log[0].format(reg=cohort_reg_method),
            ],
            "ants": [
                rules.cohort_template_reg_ants.log[0].format(reg=cohort_reg_method)
            ],
            "fireants": [],
        }[reg_presets[cohort_reg_method]["backend"]],
    shell:
        # templates built before buildinfo/ existed are archived without it
        "files=(); for f in {params.buildinfo}; do"
        "  if [ -f $f ]; then files+=($f);"
        '  else echo "Warning: $f not found, not archived" >&2; fi;'
        " done; "
        "tar -czf {output} --transform 's,.*/,,' "
        " ${{files[@]}} {params.logs}"


rule cohort_template_description:
    """template_description.json and a SPIMquant templates: config snippet."""
    input:
        files=cohort_export_files,
        build_dir=rules.build_cohort_template.output.build_dir,
    output:
        description=cohort_description,
        config_snippet=cohort_config_snippet,
    localrule: True
    params:
        name=cohort_name,
        orig_template=orig_template,
        orig_template_config=template,
        stains=cohort_stains,
        level=cohort_level,
        iterations=cohort_iterations,
        resample_mm=config["cohort_template_resample_mm"],
        subjects=[row["name"] for row in cohort_rows],
        fireants=config["cohort_template"]["fireants"],
        reg_method=cohort_reg_method,
        reg=reg_presets[cohort_reg_method],
        warps_dir=lambda wildcards, input: os.path.abspath(
            os.path.join(input.build_dir, "warps")
        ),
        orig_root=os.path.abspath(root),
        zenodo_deposition=cohort_zenodo_deposition,
        template_dir=os.path.join(cohort_root, f"tpl-{cohort_name}"),
    script:
        "../scripts/cohort_template_description.py"


def _cohort_export_file(wildcards):
    for f in cohort_upload_files:
        if os.path.basename(f) == wildcards.fname:
            return f
    raise ValueError(f"{wildcards.fname} is not a cohort template export file")


if cohort_zenodo_deposition is not None:

    rule upload_cohort_template_zenodo:
        """Upload one export file to the (draft) Zenodo deposition."""
        input:
            _cohort_export_file,
        output:
            uploaded=touch(os.path.join(cohort_work, "zenodo", "{fname}.uploaded")),
            zenodo=storage(f"zenodo://deposition/{cohort_zenodo_deposition}/{{fname}}"),
        wildcard_constraints:
            fname="[^/]+",
        localrule: True
        shell:
            "cp {input} {output.zenodo}"
