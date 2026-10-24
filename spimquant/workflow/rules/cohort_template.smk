"""
Cohort template building workflow for SPIMquant.

Builds a study-specific (cohort) SPIM template from subjects that have already
been registered to --template, then registers the original template to the new
one so its atlases can be carried over into cohort template space.

Key workflow stages:
1. Gather each subject's per-stain images in --template space (this dataset, plus
   any --cohort-template-extra-inputs SPIMquant output directories)
2. Robust intensity normalization of the subject images
3. Multivariate template building with FireANTs (all stains share one space)
4. Registration of the original template anat (moving) to the cohort template
   (fixed), since the cohort template drifts away from the original template
   space during building
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


if config.get("analysis_level") in ("cohort_template", "cohort_template_export"):
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
        bids_cohort(suffix="anat.nii.gz"),
        bids_cohort(desc="brain", suffix="mask.nii.gz") if cohort_has_mask else [],
        expand(
            bids_cohort(seg="{seg}", suffix="dseg.{ext}"),
            seg=cohort_atlas_segs,
            ext=["nii.gz", "tsv"],
        ),
        bids_cohort(from_=orig_template, to=cohort_name, suffix="xfm.nii.gz"),
        bids_cohort(from_=cohort_name, to=orig_template, suffix="xfm.nii.gz"),
        bids_cohort(
            stain=cohort_reg_stain,
            from_=orig_template,
            to=cohort_name,
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


rule cohort_template_inputs:
    """Normalize subject images and write the FireANTs image list and config.

    Each image is clipped to robust percentiles of its nonzero voxels and
    rescaled to [0,1]. The image list has one comma-separated row per subject
    (one entry per stain), which FireANTs builds as one multivariate template.
    """
    input:
        images=lambda wildcards: [img for row in cohort_rows for img in row["images"]],
    output:
        images=temp(directory(os.path.join(cohort_work, "inputs"))),
        image_list=os.path.join(cohort_work, "images.txt"),
        fireants_config=os.path.join(cohort_work, "config", "cohort.yaml"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=lambda wildcards, input: 10 + len(input.images),
    params:
        rows=cohort_rows,
        stains=cohort_stains,
        clip_percentiles=config["cohort_template"]["clip_percentiles"],
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
    """Build the multivariate cohort template with FireANTs (GPU)."""
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
        "torchrun --standalone --nproc_per_node {resources.gpu} "
        " -m fireants.scripts.template.build_template "
        " --config-dir {params.config_dir} --config-name {params.config_name} "
        " hydra.run.dir={params.hydra_dir} &> {log}"


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


rule cohort_template_affine_reg:
    """Affine registration of the original template anat to the cohort template.

    The images start out roughly aligned (subjects were registered to the
    original template), so registration is initialized with identity.
    """
    input:
        fixed=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        moving=bids_orig(suffix="anat.nii.gz"),
    output:
        xfm_ras=os.path.join(cohort_work, "xfm", "affine_ras_xfm.txt"),
        warped=os.path.join(cohort_work, "xfm", "affinewarped_anat.nii.gz"),
    log:
        os.path.join(cohort_work, "logs", "cohort_template_affine_reg.log"),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        metric=config["cohort_template"]["reg"]["metric"],
        iters=(
            "10x0x0"
            if config["sloppy"]
            else config["cohort_template"]["reg"]["affine_iters"]
        ),
    shell:
        "(greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} "
        " -a -dof 12 -ia-identity -m {params.metric} -o {output.xfm_ras} -n {params.iters} && "
        " greedy -threads {threads} -d 3 -rf {input.fixed} "
        "  -rm {input.moving} {output.warped} -r {output.xfm_ras}) &> {log}"


rule cohort_template_deform_reg:
    """Deformable registration of the original template anat to the cohort template."""
    input:
        fixed=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        moving=bids_orig(suffix="anat.nii.gz"),
        xfm_ras=rules.cohort_template_affine_reg.output.xfm_ras,
    output:
        warp=os.path.join(cohort_work, "xfm", "warp.nii.gz"),
        invwarp=os.path.join(cohort_work, "xfm", "invwarp.nii.gz"),
        warped=os.path.join(cohort_work, "xfm", "deformwarped_anat.nii.gz"),
    log:
        os.path.join(cohort_work, "logs", "cohort_template_deform_reg.log"),
    threads: 32
    resources:
        mem_mb=16000,
        runtime=5 if config["sloppy"] else 60,
    params:
        metric=config["cohort_template"]["reg"]["metric"],
        iters=(
            "10x0x0"
            if config["sloppy"]
            else config["cohort_template"]["reg"]["deform_iters"]
        ),
        sigma1=config["cohort_template"]["reg"]["sigma1"],
        sigma2=config["cohort_template"]["reg"]["sigma2"],
    shell:
        "(greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} "
        " -it {input.xfm_ras} -m {params.metric} "
        " -oinv {output.invwarp} "
        " -o {output.warp} -n {params.iters} -s {params.sigma1} {params.sigma2} && "
        " greedy -threads {threads} -d 3 -rf {input.fixed} "
        "  -rm {input.moving} {output.warped} "
        "  -r {output.warp} {input.xfm_ras}) &> {log}"


rule cohort_template_ras_to_itk:
    input:
        xfm_ras=rules.cohort_template_affine_reg.output.xfm_ras,
    output:
        xfm_itk=os.path.join(cohort_work, "xfm", "affine_itk_xfm.txt"),
    conda:
        "../envs/c3d.yaml"
    threads: 1
    resources:
        mem_mb=1500,
        runtime=15,
    shell:
        "c3d_affine_tool {input.xfm_ras} -oitk {output.xfm_itk}"


rule cohort_template_compose_warp:
    """Composite warps between the original and cohort template spaces.

    from-{orig}_to-{cohort} resamples original template images into cohort
    template space (and maps cohort-space points to the original space);
    from-{cohort}_to-{orig} is its inverse.
    """
    input:
        cohort_ref=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        orig_ref=bids_orig(suffix="anat.nii.gz"),
        xfm_itk=rules.cohort_template_ras_to_itk.output.xfm_itk,
        warp=rules.cohort_template_deform_reg.output.warp,
        invwarp=rules.cohort_template_deform_reg.output.invwarp,
    output:
        xfm_composite=bids_cohort(
            from_=orig_template, to=cohort_name, suffix="xfm.nii.gz"
        ),
        xfm_composite_inv=bids_cohort(
            from_=cohort_name, to=orig_template, suffix="xfm.nii.gz"
        ),
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


_cohort_warp_shell = (
    "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
    "antsApplyTransforms -d 3 -v -n {params.interp} "
    " -i {input.img} -o {output} -r {input.ref} -t {input.xfm}"
)


rule cohort_template_warp_anat:
    input:
        img=bids_orig(suffix="anat.nii.gz"),
        ref=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        xfm=rules.cohort_template_compose_warp.output.xfm_composite,
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
        xfm=rules.cohort_template_compose_warp.output.xfm_composite,
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
        xfm=rules.cohort_template_compose_warp.output.xfm_composite,
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

    Useful for replacing spim_templates of the original template with ones that
    are actually aligned to its atlases.
    """
    input:
        img=bids_cohort(stain="{stain}", suffix="SPIM.nii.gz"),
        ref=bids_orig(suffix="anat.nii.gz"),
        xfm=rules.cohort_template_compose_warp.output.xfm_composite_inv,
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


rule cohort_template_regqc:
    """QC report for the original-to-cohort template registration."""
    input:
        template=bids_cohort(stain=cohort_reg_stain, suffix="SPIM.nii.gz"),
        subject=bids_orig(suffix="anat.nii.gz"),
        warped_affine=rules.cohort_template_affine_reg.output.warped,
        warped_deform=rules.cohort_template_deform_reg.output.warped,
        warp=rules.cohort_template_compose_warp.output.xfm_composite,
        dseg=bids_cohort(seg=list(atlas_segs)[0], suffix="dseg.nii.gz"),
    output:
        report_html=bids_cohort(
            stain=cohort_reg_stain,
            from_=orig_template,
            to=cohort_name,
            suffix="regqc.html",
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        subject=f"tpl-{orig_template}",
        template=f"tpl-{cohort_name}",
        stain=cohort_reg_stain,
    script:
        "../scripts/reg_qc_report.py"


# --- export: prepare the cohort template for inclusion in SPIMquant ---

# every file a SPIMquant templates: entry for the cohort template refers to
cohort_export_files = (
    [bids_cohort(suffix="anat.nii.gz")]
    + ([bids_cohort(desc="brain", suffix="mask.nii.gz")] if cohort_has_mask else [])
    + [bids_cohort(stain=s, suffix="SPIM.nii.gz") for s in cohort_stains]
    + [
        bids_cohort(seg=seg, suffix=f"dseg.{ext}")
        for seg in cohort_atlas_segs
        for ext in ["nii.gz", "tsv"]
    ]
    + [
        bids_cohort(from_=orig_template, to=cohort_name, suffix="xfm.nii.gz"),
        bids_cohort(from_=cohort_name, to=orig_template, suffix="xfm.nii.gz"),
    ]
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
        image_list=rules.cohort_template_inputs.output.image_list,
        fireants_config=rules.cohort_template_inputs.output.fireants_config,
        # logs cannot be inputs; these outputs guarantee the logs were written
        build_dir=rules.build_cohort_template.output.build_dir,
        warp=rules.cohort_template_deform_reg.output.warp,
    output:
        cohort_buildinfo,
    localrule: True
    params:
        logs=[
            rules.build_cohort_template.log[0],
            rules.cohort_template_affine_reg.log[0],
            rules.cohort_template_deform_reg.log[0],
        ],
    shell:
        "tar -czf {output} --transform 's,.*/,,' "
        " {input.image_list} {input.fireants_config} {params.logs}"


rule cohort_template_description:
    """template_description.json and a SPIMquant templates: config snippet."""
    input:
        files=cohort_export_files,
        image_list=rules.cohort_template_inputs.output.image_list,
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
        subjects=[row["name"] for row in cohort_rows],
        fireants=config["cohort_template"]["fireants"],
        reg=config["cohort_template"]["reg"],
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
