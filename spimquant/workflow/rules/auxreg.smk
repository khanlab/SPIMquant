"""
Registration of additional scans (--aux_spim) to their reference scan.

An additional scan (e.g. a 12x scan alongside the 1x scan of a sample) is not
registered to the template itself: it is registered to the reference scan of the
same sample, and its transforms to and from the template chain that registration
with the reference scan's. They are written where the template registration
methods write theirs, with desc-viaref, so select_subject_to_template_warp picks
them up (see subject_reg_method) and everything downstream (atlas labels, template
space outputs and quantification) works for additional scans as it does for any
scan.

The additional scan is the fixed image and the reference scan the moving image,
so the similarity metric is computed over the additional scan's field of view,
which may cover only part of the sample. All of --aux_reg_stains are registered
together (multichannel), on the scans as imported (without bias correction or
masking, since the additional scan's brain mask is derived from this
registration, which is why these rules are in the subj_mask group).

Transforms (greedy convention: a transform maps points in the fixed image to the
moving image, and -r applies the transforms listed in order):

    G = -r [warp] affine          additional scan -> reference scan (from-ref_to-subject)
    G^-1 = -r affine,-1 [invwarp] reference scan -> additional scan

    from-subject_to-{template}: -r <reference scan's from-subject_to-{template}> G^-1
    from-{template}_to-subject: -r G <reference scan's from-{template}_to-subject>
"""

aux_deform = config["aux_reg_method"] == "deform"


def _aux_reg_image(stain):
    return bids(
        root=root,
        datatype="micr",
        stain=stain,
        level=config["registration_level"],
        suffix="SPIM.nii.gz",
        **inputs["spim"].wildcards,
    )


def _aux_xfm(suffix, from_="ref", to="subject", **kwargs):
    return bids(
        root=root,
        datatype="xfm",
        from_=from_,
        to=to,
        suffix=suffix,
        **kwargs,
        **inputs["spim"].wildcards,
    )


def _aux_ref_images(wildcards):
    return [
        _aux_reg_image(stain).format(**ref_scan_wildcards(wildcards))
        for stain in aux_reg_stains
    ]


def _aux_to_ref(input):
    """greedy -r transforms from the additional scan to its reference scan (G)."""
    return f"{input.warp} {input.xfm_ras}" if aux_deform else f"{input.xfm_ras}"


def _aux_reg_cmd(wildcards, input, output, threads):
    """greedy commands for aux_reg_to_ref: rigid, then affine, then deformable."""
    opts = config["aux_reg"]
    greedy = f"greedy -threads {threads} -d 3 " + " ".join(
        f"-i {f} {m}" for f, m in zip(input.fixed, input.moving)
    )
    metric = f"-m {opts['metric']}"
    cmds = [
        f"{greedy} -a -dof 6 -ia {input.init} {metric} -n {opts['rigid_iters']} "
        f"-o {output.xfm_ras}"
    ]
    if config["aux_reg_method"] != "rigid":
        cmds.append(
            f"{greedy} -a -dof 12 -ia {output.xfm_ras} {metric} "
            f"-n {opts['affine_iters']} -o {output.xfm_ras}"
        )
    if aux_deform:
        cmds.append(
            f"{greedy} -it {output.xfm_ras} {metric} -n {opts['deform_iters']} "
            f"-s {opts['sigma1']} {opts['sigma2']} "
            f"-o {output.warp} -oinv {output.invwarp}"
        )
    cmds.append(
        f"greedy -threads {threads} -d 3 -rf {input.fixed[0]} "
        f"-rm {input.moving[0]} {output.warped} -r {_aux_to_ref(output)}"
    )
    return " && ".join(cmds)


rule aux_reg_init:
    """Initial translation of an additional scan within its reference scan.

    See --aux_reg_init: by default, a normalized cross-correlation search for the
    additional scan's position within the reference scan, as the physical
    coordinates of two scans need not relate to each other.
    """
    input:
        fixed=expand(
            _aux_reg_image("{stain}"), stain=aux_reg_stains, allow_missing=True
        ),
        moving=_aux_ref_images,
    output:
        xfm_ras=temp(_aux_xfm("xfm.txt", desc="init", type_="ras")),
        # the translation and match score, for aux_reg_qc_report
        info=temp(_aux_xfm("xfm.json", desc="init")),
    log:
        bids(
            root="logs",
            datatype="aux_reg_init",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_mask"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        method=config["aux_reg_init"],
        min_res_um=config["aux_reg"]["match_min_res_um"],
    script:
        "../scripts/aux_reg_init.py"


rule aux_reg_to_ref:
    """Register an additional scan to its reference scan with greedy (--aux_reg_method).

    Rigid, then affine unless --aux_reg_method rigid, then deformable with
    --aux_reg_method deform, all multichannel over --aux_reg_stains. Also writes
    the reference scan resampled to the additional scan, to check the alignment.
    """
    input:
        fixed=expand(
            _aux_reg_image("{stain}"), stain=aux_reg_stains, allow_missing=True
        ),
        moving=_aux_ref_images,
        init=rules.aux_reg_init.output.xfm_ras,
    output:
        **(
            {
                "warp": temp(_aux_xfm("warp.nii.gz")),
                "invwarp": temp(_aux_xfm("warp.nii.gz", from_="subject", to="ref")),
            }
            if aux_deform
            else {}
        ),
        xfm_ras=_aux_xfm("xfm.txt", desc="affine", type_="ras"),
        warped=bids(
            root=root,
            datatype="xfm",
            space="subject",
            from_="ref",
            stain=stain_for_reg,
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
    log:
        bids(
            root="logs",
            datatype="aux_reg_to_ref",
            suffix="log.txt",
            **inputs["spim"].wildcards,
        ),
    group:
        "subj_mask"
    threads: 32
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        cmd=_aux_reg_cmd,
    shell:
        "({params.cmd}) &> {log}"


rule compose_aux_to_template_via_ref:
    """Additional scan <-> template warps, through the reference scan (desc-viaref).

    Not in the subj_reg group: it depends on the reference scan's subj_reg jobs
    through select_subject_to_template_warp (a local rule), which would make the
    additional scan's group depend on itself.
    """
    input:
        **(
            {
                "warp": rules.aux_reg_to_ref.output.warp,
                "invwarp": rules.aux_reg_to_ref.output.invwarp,
            }
            if aux_deform
            else {}
        ),
        ref=rules.import_template_anat.output.anat,
        subject=_aux_reg_image(stain_for_reg),
        xfm_ras=rules.aux_reg_to_ref.output.xfm_ras,
        ref_to_template=ref_scan_input(
            rules.select_subject_to_template_warp.output.xfm_composite
        ),
        template_to_ref=ref_scan_input(
            rules.select_subject_to_template_warp.output.xfm_composite_inv
        ),
    output:
        xfm_composite=_reg_xfm("xfm.nii.gz", desc="viaref"),
        xfm_composite_inv=_reg_xfm(
            "xfm.nii.gz", from_="{template}", to="subject", desc="viaref"
        ),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=15,
    params:
        to_ref=lambda wildcards, input: _aux_to_ref(input),
        from_ref=lambda wildcards, input: (
            f"{input.xfm_ras},-1 {input.invwarp}"
            if aux_deform
            else f"{input.xfm_ras},-1"
        ),
    shell:
        "greedy -threads {threads} -d 3 -rf {input.ref} -rc {output.xfm_composite} "
        " -r {input.ref_to_template} {params.from_ref} && "
        "greedy -threads {threads} -d 3 -rf {input.subject} -rc {output.xfm_composite_inv} "
        " -r {params.to_ref} {input.template_to_ref}"


rule warp_ref_mask_to_aux:
    """Reference scan's brain mask in an additional scan, as the prior for its mask.

    Takes the place of the template mask prior (affine_transform_template_mask_to_subject)
    for additional scans, whose field of view may not hold enough of the sample
    for an affine registration to the template.
    """
    input:
        **({"warp": rules.aux_reg_to_ref.output.warp} if aux_deform else {}),
        ref=_aux_reg_image(stain_for_reg),
        mask=ref_scan_input(
            bids(
                root=root,
                datatype="micr",
                stain=stain_for_reg,
                level=config["registration_level"],
                desc="brain",
                suffix="mask.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
        xfm_ras=rules.aux_reg_to_ref.output.xfm_ras,
    output:
        mask=temp(
            bids(
                root=root,
                datatype="micr",
                desc="viaref",
                from_="ref",
                suffix="mask.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj_mask"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=10,
    params:
        to_ref=lambda wildcards, input: _aux_to_ref(input),
    shell:
        "greedy -threads {threads} -d 3 -rf {input.ref} -ri NN "
        " -rm {input.mask} {output.mask} -r {params.to_ref}"


def _scan_name(key):
    """BIDS name of a scan, from its entity tuple (_spim_key_wildcards)."""
    return os.path.basename(inputs["spim"].path).format(
        **dict(zip(_spim_key_wildcards, key))
    )


def _ref_to_aux(input):
    """greedy -r transforms from the reference scan to its additional scan (G^-1)."""
    return (
        f"{input.xfm_ras},-1 {input.invwarp}" if aux_deform else f"{input.xfm_ras},-1"
    )


rule aux_coverage_tsv:
    """Fraction of each atlas region that an additional scan's field of view covers.

    The field of view (the additional scan's nonzero voxels) is resampled into its
    reference scan, and averaged over each region of the reference scan's atlas
    labels, which cover the whole sample: the coverage column of the additional
    scan's mergedsegstats.tsv.
    """
    input:
        **({"invwarp": rules.aux_reg_to_ref.output.invwarp} if aux_deform else {}),
        fov=_aux_reg_image(stain_for_reg),
        xfm_ras=rules.aux_reg_to_ref.output.xfm_ras,
        dseg=ref_scan_input(
            bids(
                root=root,
                datatype="parc",
                seg="{seg}",
                level=config["registration_level"],
                from_="{template}",
                suffix="dseg.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
        label_tsv=bids(root=root, template="{template}", seg="{seg}", suffix="dseg.tsv"),
    output:
        tsv=temp(
            bids(
                root=root,
                datatype="tabular",
                seg="{seg}",
                from_="{template}",
                suffix="coveragestats.tsv",
                **inputs["spim"].wildcards,
            )
        ),
    group:
        "subj_quant"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=10,
    params:
        to_aux=lambda wildcards, input: _ref_to_aux(input),
    script:
        "../scripts/aux_coverage.py"


# the reference scan resampled into an additional scan, with the initial
# transform (init) or the registration (final), for aux_reg_qc_report
_aux_qc_warped = bids(
    root=root,
    datatype="xfm",
    from_="ref",
    stain="{stain}",
    desc="{stage}",
    suffix="auxqc.nii.gz",
    **inputs["spim"].wildcards,
)


rule warp_ref_to_aux_for_report:
    """Resample a stain of the reference scan into an additional scan, for the report."""
    input:
        **({"warp": rules.aux_reg_to_ref.output.warp} if aux_deform else {}),
        fixed=_aux_reg_image("{stain}"),
        moving=ref_scan_input(_aux_reg_image("{stain}")),
        init=rules.aux_reg_init.output.xfm_ras,
        xfm_ras=rules.aux_reg_to_ref.output.xfm_ras,
    output:
        temp(_aux_qc_warped),
    wildcard_constraints:
        stage="init|final",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=10,
    params:
        xfms=lambda wildcards, input: (
            input.init if wildcards.stage == "init" else _aux_to_ref(input)
        ),
    shell:
        "greedy -threads {threads} -d 3 -rf {input.fixed} "
        " -rm {input.moving} {output} -r {params.xfms}"


rule aux_reg_qc_report:
    """QC report of an additional scan's registration to its reference scan.

    Shows the reference scan resampled into the additional scan with the initial
    translation (--aux_reg_init) and with the registration (--aux_reg_method), for
    each of --aux_reg_stains, with the atlas labels the additional scan gets
    through its reference scan.
    """
    input:
        **({"warp": rules.aux_reg_to_ref.output.warp} if aux_deform else {}),
        fixed=expand(
            _aux_reg_image("{stain}"), stain=aux_reg_stains, allow_missing=True
        ),
        warped_init=expand(
            _aux_qc_warped, stain=aux_reg_stains, stage="init", allow_missing=True
        ),
        warped_final=expand(
            _aux_qc_warped, stain=aux_reg_stains, stage="final", allow_missing=True
        ),
        dseg=bids(
            root=root,
            datatype="parc",
            seg=list(atlas_segs)[0],
            level=config["registration_level"],
            from_=config["template"],
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        init_info=rules.aux_reg_init.output.info,
    output:
        report_html=_aux_xfm("regqc.html"),
    threads: 1
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        scan=lambda wildcards: _scan_name(_scan_key(wildcards)),
        ref=lambda wildcards: _scan_name(aux_ref_keys[_scan_key(wildcards)]),
        stains=aux_reg_stains,
        init=config["aux_reg_init"],
        method=config["aux_reg_method"],
        # those of the stages --aux_reg_method runs
        parameters=" ".join(
            f"{k}={config['aux_reg'][k]}"
            for k in ["metric", "rigid_iters"]
            + (["affine_iters"] if config["aux_reg_method"] != "rigid" else [])
            + (["deform_iters", "sigma1", "sigma2"] if aux_deform else [])
        ),
    script:
        "../scripts/aux_reg_qc_report.py"
