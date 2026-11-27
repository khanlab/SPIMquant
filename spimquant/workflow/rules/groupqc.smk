"""
Group QC report: registration and segmentation QC of all subjects in one page.

The participant level writes QC figures and reports for each subject; this
gathers them at the group level. Each subject's stains and field fraction maps
are already resampled to template space, so the report shows the same template
slices for every subject: one subject at a time with its metrics against the
group, a montage of all subjects at one slice, and an animated flip through
them. Subject-level figures (registration report, segmentation slices, ROI
montages, ...) are linked from the subject view, or embedded with
--groupqc_embed_figures.

groupqc_subject renders the slices and metrics of one subject (light, run per
subject); groupqc_report adds what needs the whole group (outlier flags,
correlation with the group mean, a shared color scale) and writes the page and
a table of the metrics.
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

groupqc_root = os.path.join(root, "group", "qc")


def _groupqc_outline_seg():
    if config["groupqc"].get("outline_seg"):
        return config["groupqc"]["outline_seg"]
    for seg in ("coarse", "roi22"):
        if seg in template.get("atlases", {}):
            return seg
    return list(atlas_segs)[0] if atlas_segs else None


groupqc_outline_seg = _groupqc_outline_seg()

# the registration stain first: it is the image the registration is checked on
groupqc_stains = [stain_for_reg] + (
    [s for s in stains if s != stain_for_reg] if do_seg else []
)

# field fraction maps in template space: each segmented stain's primary method,
# and the vessels
groupqc_seg_layers = [
    {"stain": p.stain, "desc": p.seg_method, "datatype": "seg"} for p in seg_plan
] + (
    [
        {
            "stain": stain_for_vessels,
            "desc": config["vessel_seg_method"],
            "datatype": "vessels",
        }
    ]
    if do_vessels
    else []
)


def _groupqc_template_ref():
    """Template image the registration stain is compared with."""
    spim_templates = template.get("spim_templates") or {}
    suffix = (
        stain_for_reg
        if stain_for_reg in spim_templates and config["template"] == reg_template
        else "anat"
    )
    return bids(root=root, template=config["template"], suffix=f"{suffix}.nii.gz")


def _groupqc_subject_prefix(datatype):
    """Path prefix of a subject's files in a datatype folder, for globbing."""

    def prefix(wildcards):
        path = bids(
            root=root, datatype=datatype, suffix="x", **inputs["spim"].wildcards
        ).format(**wildcards)
        return path[: -len("_x")]

    return prefix


rule groupqc_subject:
    """Template-space slices and QC metrics of one subject, for the group report."""
    input:
        template_anat=_groupqc_template_ref(),
        template_mask=bids(
            root=root, template=config["template"], desc="brain", suffix="mask.nii.gz"
        ),
        warped=[
            bids(
                root=root,
                datatype="micr",
                stain=stain,
                level=config["registration_level"],
                space=config["template"],
                suffix="SPIM.nii.gz",
                **inputs["spim"].wildcards,
            )
            for stain in groupqc_stains
        ],
        fieldfrac=[
            bids(
                root=root,
                datatype=layer["datatype"],
                stain=layer["stain"],
                level=config["registration_level"],
                desc=layer["desc"],
                space=config["template"],
                suffix="fieldfrac.nii.gz",
                **inputs["spim"].wildcards,
            )
            for layer in groupqc_seg_layers
        ],
        xfm=(
            bids(
                root=root,
                datatype="xfm",
                from_="subject",
                to=config["template"],
                suffix="xfm.nii.gz",
                **inputs["spim"].wildcards,
            )
            if config["groupqc"]["jacobian"]
            else []
        ),
        segstats=(
            bids(
                root=root,
                datatype="tabular",
                seg=list(atlas_segs)[0],
                from_=config["template"],
                suffix="mergedsegstats.tsv",
                **inputs["spim"].wildcards,
            )
            if do_seg
            else []
        ),
    output:
        json=bids(
            root=root, datatype="qc", suffix="groupqc.json", **inputs["spim"].wildcards
        ),
    group:
        "subj_groupqc"
    threads: 1
    resources:
        mem_mb=16000,
        runtime=20,
    params:
        label=lambda wildcards: "_".join(
            f"{k}-{v}" if k != "subject" else f"sub-{v}"
            for k, v in wildcards.items()
            if k == "subject" or len({str(x) for x in inputs["spim"].zip_lists[k]}) > 1
        ),
        entities=lambda wildcards: dict(wildcards),
        stains=groupqc_stains,
        seg_layers=groupqc_seg_layers,
        n_slices=config["groupqc"]["slices"],
        quality=config["groupqc"]["image_quality"],
        figure_prefix=_groupqc_subject_prefix("qc"),
        regqc=lambda wildcards: [
            _groupqc_subject_prefix("xfm")(wildcards) + "_*regqc.html"
        ],
    script:
        "../scripts/groupqc_subject.py"


def _groupqc_participants():
    path = os.path.join(str(config["bids_dir"]), "participants.tsv")
    return path if os.path.exists(path) else []


def _groupqc_version():
    try:
        return _pkg_version("spimquant")
    except PackageNotFoundError:
        return ""


rule groupqc_report:
    """One self-contained HTML page with the QC of every subject, and a metrics table."""
    input:
        subjects=ref_spim.expand(rules.groupqc_subject.output.json),
        template_anat=_groupqc_template_ref(),
        template_mask=rules.groupqc_subject.input.template_mask,
        outline_dseg=(
            bids(
                root=root,
                template=config["template"],
                seg=groupqc_outline_seg,
                suffix="dseg.nii.gz",
            )
            if groupqc_outline_seg
            else []
        ),
        participants=_groupqc_participants(),
    output:
        html=os.path.join(groupqc_root, "groupqc.html"),
        tsv=os.path.join(groupqc_root, "groupqc.tsv"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        title=os.path.basename(os.path.normpath(str(config["bids_dir"]))),
        template=config["template"],
        group_by=config.get("groupqc_group_by")
        or config.get("group_stats_pairwise")
        or [],
        outlier_z=config["groupqc"]["outlier_z"],
        embed_figures=config.get("groupqc_embed_figures", False),
        reg_method=config["registration_method"],
        version=_groupqc_version(),
        quality=config["groupqc"]["image_quality"],
    script:
        "../scripts/groupqc_report.py"


rule all_groupqc:
    input:
        rules.groupqc_report.output,
