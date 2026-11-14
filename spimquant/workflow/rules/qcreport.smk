"""
Dataset QC report on the raw data (the qc analysis level).

A light alternative to the participant level, to check a batch of scans before
committing to segmentation: per-scan intensity, coverage and sharpness metrics,
overview slices and full-resolution crops of every stain, gathered in one
self-contained HTML report with selectors for scan and stain.

Unlike the rest of the workflow, scans do not need to have the same stains:
each scan is reported on with the stains it has (qc_scan_stains), so datasets
that pilot different stainings can be compared. For this reason these rules do
not use the workflow's registration: each scan is registered to the template
with a quick greedy registration of its own first stain from --stains_for_reg,
on the raw image, only to place the crops in the same atlas regions in every
scan. Scans without such a stain (or all, with --qc_no_reg) get crops spread
through the tissue instead.
"""

qc_template = config["templates"][config["template"]]
qc_reg_preset = reg_presets[config["qc"]["reg_preset"]]

if config["qc_atlas_seg"] is not None:
    if config["qc_atlas_seg"] not in qc_template.get("atlases", {}):
        raise ValueError(
            f"--qc_atlas_seg {config['qc_atlas_seg']} was not found in the template "
            f"{config['template']}"
        )
    qc_atlas_seg = config["qc_atlas_seg"]
elif "coarse" in qc_template.get("atlases", {}):
    qc_atlas_seg = "coarse"
else:
    qc_atlas_seg = next(iter(qc_template.get("atlases", {})), None)


def qc_stains(wildcards):
    """Stains of a scan, in channel order."""
    return qc_scan_stains[str(spim_input(wildcards))]


def qc_reg_stain(wildcards):
    """The stain a scan is registered with, or None if it is not registered."""
    if config["qc_no_reg"] or qc_atlas_seg is None:
        return None
    scan_stains = qc_stains(wildcards)
    return next((s for s in config["stains_for_reg"] if s in scan_stains), None)


def _qc_require_reg_stain(wildcards):
    stain = qc_reg_stain(wildcards)
    if stain is None:
        raise ValueError(f"No stain to register this scan with: {dict(wildcards)}")
    return stain


def qc_reg_fixed(wildcards):
    """Template image of the scan's registration stain, or else its anat."""
    stain = _qc_require_reg_stain(wildcards)
    suffix = stain if stain in (qc_template.get("spim_templates") or {}) else "anat"
    return get_template_path(
        root, config["template"], config.get("template_crop"), suffix
    )


def qc_reg_moving(wildcards):
    return bids(
        root=root,
        datatype="micr",
        stain=_qc_require_reg_stain(wildcards),
        level=config["registration_level"],
        suffix="SPIM.nii.gz",
        **inputs["spim"].wildcards,
    ).format(**wildcards)


def _qc_out(suffix, **kwargs):
    return bids(
        root=root, datatype="qc", suffix=suffix, **kwargs, **inputs["spim"].wildcards
    )


rule qc_register:
    """Quick registration of a scan to the template, for the qc report.

    Affine then deformable greedy registration of the raw downsampled image
    (no bias field correction or masking), and the atlas labels in the space of
    the scan.
    """
    input:
        fixed=qc_reg_fixed,
        moving=qc_reg_moving,
        dseg=bids(
            root=root,
            template=config["template"],
            seg=qc_atlas_seg or "none",
            suffix="dseg.nii.gz",
        ),
    output:
        xfm_ras=temp(_qc_out("xfm.txt", desc="qcreg", from_="subject", to="template")),
        warp=temp(_qc_out("warp.nii.gz", desc="qcreg", from_="subject", to="template")),
        invwarp=temp(
            _qc_out("warp.nii.gz", desc="qcreg", from_="template", to="subject")
        ),
        warped=temp(_qc_out("SPIM.nii.gz", desc="qcreg", space=config["template"])),
        dseg=_qc_out("dseg.nii.gz", desc="qcreg", from_=config["template"]),
    log:
        bids(
            root="logs",
            datatype="qc_register",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_qc"
    threads: 32
    resources:
        mem_mb=per_attempt(16000),
        runtime=5 if config["sloppy"] else 30,
    params:
        metric=qc_reg_preset["metric"],
        affine_iters=reg_iters(qc_reg_preset, "affine_iters"),
        deform_iters=reg_iters(qc_reg_preset, "deform_iters"),
        sigma1=qc_reg_preset["sigma1"],
        sigma2=qc_reg_preset["sigma2"],
    shell:
        "(greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} "
        "  -a -dof 12 -ia-image-centers -m {params.metric} -o {output.xfm_ras} "
        "  -n {params.affine_iters} && "
        " greedy -threads {threads} -d 3 -i {input.fixed} {input.moving} "
        "  -it {output.xfm_ras} -m {params.metric} -o {output.warp} "
        "  -oinv {output.invwarp} -n {params.deform_iters} "
        "  -s {params.sigma1} {params.sigma2} && "
        " greedy -threads {threads} -d 3 -rf {input.fixed} "
        "  -rm {input.moving} {output.warped} -r {output.warp} {output.xfm_ras} && "
        " greedy -threads {threads} -d 3 -rf {input.moving} -ri NN "
        "  -rm {input.dseg} {output.dseg} -r {output.xfm_ras},-1 {output.invwarp}"
        ") &> {log}"


def _qc_scan_reg_inputs(wildcards):
    if qc_reg_stain(wildcards) is None:
        return {}
    return {
        "dseg": rules.qc_register.output.dseg.format(**wildcards),
        "warped": rules.qc_register.output.warped.format(**wildcards),
        "fixed": qc_reg_fixed(wildcards),
        "template_dseg": rules.qc_register.input.dseg,
        "label_tsv": bids(
            root=root,
            template=config["template"],
            seg=qc_atlas_seg,
            suffix="dseg.tsv",
        ),
    }


rule qc_scan:
    """QC metrics and images of one scan, for the qc report.

    For each stain of the scan: intensity, contrast and depth-uniformity
    metrics and overview slices from the downsampled image, and
    full-resolution crops, with their sharpness and saturation, in atlas
    regions (or spread through the tissue, if the scan is not registered).
    Written as one JSON with the images embedded.
    """
    input:
        unpack(_qc_scan_reg_inputs),
        spim=spim_input,
        spim_json=spim_json_input,
    output:
        json=_qc_out("qcscan.json"),
    group:
        "subj_qc"
    threads: 8
    resources:
        mem_mb=scaled_mem_mb(8000, 16000, config["registration_level"]),
        runtime=30,
    params:
        stains=qc_stains,
        reg_stain=qc_reg_stain,
        entities=lambda wildcards: dict(wildcards),
        template=config["template"],
        atlas_seg=qc_atlas_seg,
        level=config["registration_level"],
        crop_level=config["qc_crop_level"],
        crop_size=config["qc_crop_size"],
        n_crops=config["qc_n_crops"],
        crop_labels=config["qc_crop_labels"],
        overview_slices=config["qc"]["overview_slices"],
        overview_size=config["qc"]["overview_size"],
        image_quality=config["qc"]["image_quality"],
        itemsize=lambda wildcards, input: get_spim_info(input.spim, 0)["itemsize"],
        zarrnii_kwargs=zarrnii_in_kwargs,
    script:
        "../scripts/qc_scan.py"


_qc_participants_tsv = os.path.join(config["bids_dir"], "participants.tsv")


rule qc_dataset_report:
    """One self-contained HTML report on all scans, and a table of their metrics."""
    input:
        scans=inputs["spim"].expand(rules.qc_scan.output.json),
        participants=(
            _qc_participants_tsv if os.path.exists(_qc_participants_tsv) else []
        ),
    output:
        html=os.path.join(root, "qc", "qcreport.html"),
        tsv=os.path.join(root, "qc", "qcmetrics.tsv"),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    params:
        title=os.path.basename(os.path.normpath(config["bids_dir"])),
        stain_aliases=config["qc"]["stain_aliases"],
        group_by=config["qc_group_by"],
    script:
        "../scripts/qc_dataset_report.py"


rule all_qcreport:
    input:
        rules.qc_dataset_report.output,
