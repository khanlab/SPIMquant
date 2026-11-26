"""
Landmark evaluation of registration (--landmark_analysis).

Registers every SPIM scan to --template with every registration method
(--landmark_methods, default: every single-channel registration.presets entry)
using every stain as the registration channel (--landmark_stains, default: all
stains of the scans, ignoring the --stains_for_reg priority), and with
--register_to_mri registers the MRI to SPIM once per stain. Then computes the
target registration error (TRE) of manual landmarks for every transform:

  spim_to_template   SPIM landmarks -> template landmarks
  template_to_spim   template landmarks -> SPIM landmarks
  mri_to_spim        MRI landmarks -> SPIM landmarks          (per stain)
  spim_to_mri        SPIM landmarks -> MRI landmarks          (per stain)
  mri_to_template    MRI -> SPIM -> template, chained          (per stain x method)

Points move opposite to images, so e.g. SPIM points go to the template with the
template->subject transform. SPIM->template is reported twice: with the shipped
inverse composite (inversion=shipped) and by fixed-point inversion of the
forward composite (inversion=fixedpoint).

The registration rules reuse the regular ones (use rule ... with:), with the
registration stain as the {regstain} wildcard and outputs under landmarks/.
The brain mask (and so the masked N4 image of every stain) is the one of the
default registration stain, so stains differ only in the registration channel.
Scans may have different stains: each is registered with the stains it has, and
the brain mask stain is one that every scan has.

Besides the TRE tables, every warped landmark set is written as a 3D Slicer
markups file (landmark_markups), so where the points land can be looked at.
"""

import glob as _glob
import sys as _sys
from types import SimpleNamespace

_sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
import landmark_files as LF

lmk_root = os.path.join(root, "landmarks")

lmk_methods = config.get("landmark_methods") or [
    n for n, p in reg_presets.items() if not p.get("multichannel")
]
for _m in lmk_methods:
    if _m not in reg_presets:
        raise ValueError(
            f"--landmark_methods: unknown method {_m}, choose from {sorted(reg_presets)}"
        )
    if reg_presets[_m].get("multichannel"):
        raise ValueError(
            f"--landmark_methods: {_m} is multichannel, not supported here"
        )

lmk_stains = config.get("landmark_stains") or list(stains)
for _s in lmk_stains:
    if _s not in stains:
        raise ValueError(f"--landmark_stains: {_s} not among the scan stains {stains}")


def lmk_stains_of(wildcards):
    """The registration stains of one scan: --landmark_stains that the scan has."""
    if not lmk_scan_stains:
        return lmk_stains
    w = dict(wildcards.items()) if hasattr(wildcards, "items") else vars(wildcards)
    zl = spim_raw.zip_lists
    keys = list(zl)
    for vals in zip(*[zl[k] for k in keys]):
        row = dict(zip(keys, vals))
        if all(str(w[k]) == str(row[k]) for k in keys if k in w):
            have = lmk_scan_stains.get(spim_raw.path.format(**row), lmk_stains)
            out = [s for s in lmk_stains if s in have]
            if not out:
                raise ValueError(
                    f"--landmark_analysis: scan {w} has none of the stains {lmk_stains}"
                )
            return out
    return lmk_stains


lmk_do_mri = bool(config["register_to_mri"] and config.get("landmarks_mri"))
if config.get("landmarks_mri") and not config["register_to_mri"]:
    raise ValueError("--landmarks_mri needs --register_to_mri")
for _k in ("landmarks_spim", "landmarks_template"):
    if not config.get(_k):
        raise ValueError(f"--landmark_analysis needs --{_k}")

_lmk_gpu_methods = [
    m for m in lmk_methods if reg_presets[m].get("backend") == "fireants"
]
if _lmk_gpu_methods and not gpu_available:
    # gpu_available (Snakefile): --gpu, or by default whether torch has CUDA
    print(
        f"Landmark analysis: NOTE {', '.join(_lmk_gpu_methods)} need(s) a GPU (the gpu pixi "
        "environment and a GPU for the jobs), and this environment has none (--gpu). If the "
        "jobs will not get one, leave it out with: --landmark_methods "
        + " ".join(m for m in lmk_methods if m not in _lmk_gpu_methods)
    )

print(
    f"Landmark analysis: methods {lmk_methods}, registration stains {lmk_stains}"
    + (", with MRI" if lmk_do_mri else "")
)
if len({tuple(v) for v in lmk_scan_stains.values()}) > 1:
    print(
        f"Landmark analysis: scans have different stains, each uses its own: {lmk_scan_stains}"
    )


wildcard_constraints:
    regstain="|".join(lmk_stains),


# ---------------------------------------------------------------- landmark files


def _lmk_find(base, wildcards, datatype, require=None):
    """The landmarks file of a scan in a manual_landmarks derivative, or None."""
    subj = f"sub-{wildcards.subject}"
    hits = sorted(
        _glob.glob(
            os.path.join(base, subj, "**", datatype, f"{subj}_*landmarks.json"),
            recursive=True,
        )
    )
    if require:
        hits = [h for h in hits if require in os.path.basename(h)]
    # a file naming another session/sample/acquisition belongs to another scan
    return LF.pick_landmark_file(hits, LF.scan_entities(wildcards))


def _mri_sessions(subject):
    try:
        return sorted(
            set(inputs["mri"].filter(subject=subject).entities.get("session", []))
        )
    except Exception:
        return []


def _lmk_mri_entities(path):
    """{MRI wildcard: value} named in a landmark file name, for wildcards the MRI input has."""
    return LF.mri_entities(path, inputs["mri"].wildcards)


def lmk_spim_file(wildcards):
    return _lmk_find(config["landmarks_spim"], wildcards, "micr", require="from-SPIM")


def lmk_mri_file(wildcards):
    if not lmk_do_mri:
        return None
    base, subj = config["landmarks_mri"], f"sub-{wildcards.subject}"
    hits = sorted(
        _glob.glob(
            os.path.join(
                base, subj, "**", "anat", f"{subj}_*from-{mri_suffix}_landmarks.json"
            ),
            recursive=True,
        )
    )
    sessions = _mri_sessions(wildcards.subject)
    if len(hits) > 1 and sessions:
        hits = [
            h
            for h in hits
            if any(f"_ses-{s}_" in os.path.basename(h) for s in sessions)
        ]
    if len(hits) > 1:
        raise ValueError(
            f"several MRI landmark files for sub-{wildcards.subject}: {hits}; "
            "use --filter-mri session=<label> to choose one"
        )
    return hits[0] if hits else None


# ---------------------------------------------------------------- paths


def _lmk_moving(level=None):
    return bids(
        root=root,
        datatype="micr",
        stain="{regstain}",
        level=level if level is not None else config["registration_level"],
        desc=config["templatereg"]["desc"],
        suffix="SPIM.nii.gz",
        **inputs["spim"].wildcards,
    )


def _lmk_xfm(suffix, from_="subject", to="{template}", desc="{reg}", **kwargs):
    ents = dict(stain="{regstain}", from_=from_, to=to, desc=desc, **kwargs)
    return bids(
        root=lmk_root,
        datatype="xfm",
        suffix=suffix,
        **{k: v for k, v in ents.items() if v is not None},
        **inputs["spim"].wildcards,
    )


def _lmk_warped(desc, space="{template}", suffix="SPIM.nii.gz", **kwargs):
    return bids(
        root=lmk_root,
        datatype="xfm",
        space=space,
        stain="{regstain}",
        desc=desc,
        suffix=suffix,
        **kwargs,
        **inputs["spim"].wildcards,
    )


def _lmk_log(name, **kwargs):
    return bids(
        root="logs",
        datatype=name,
        stain="{regstain}",
        suffix="log.txt",
        **kwargs,
        **inputs["spim"].wildcards,
    )


_fwd = _lmk_xfm("xfm.nii.gz")
lmk_xfm_fwd = temp(_fwd) if not config.get("landmark_keep_warps") else _fwd
lmk_xfm_inv = _lmk_xfm("xfm.nii.gz", from_="{template}", to="subject")


# ---------------------------------------------------------------- SPIM <-> template


use rule affine_reg as lmk_affine_reg with:
    input:
        fixed=get_template_for_reg,
        moving=_lmk_moving(),
        extra_fixed=[],
        extra_moving=[],
    output:
        xfm_ras=temp(_lmk_xfm("xfm.txt", desc="{reg}affine", type_="ras")),
        warped=temp(_lmk_warped("{reg}affinewarped")),
    log:
        _lmk_log("lmk_affine_reg", space="{template}", desc="{reg}"),
    group:
        "lmk_reg"
    params:
        affine=lambda wildcards: reg_presets[wildcards.reg]["affine"],
        init="-ia-image-centers",
        metric=lambda wildcards: reg_presets[wildcards.reg]["metric"],
        iters=lambda wildcards: reg_iters(reg_presets[wildcards.reg], "affine_iters"),
        identity=reg_identity_ras,
        extra_pairs="",


use rule convert_ras_to_itk as lmk_convert_ras_to_itk with:
    input:
        xfm_ras=rules.lmk_affine_reg.output.xfm_ras,
    output:
        xfm_itk=temp(_lmk_xfm("xfm.txt", desc="{reg}affine", type_="itk")),
    group:
        "lmk_reg"


use rule deform_reg as lmk_deform_reg with:
    input:
        fixed=get_template_for_reg,
        moving=_lmk_moving(),
        extra_fixed=[],
        extra_moving=[],
        xfm_ras=rules.lmk_affine_reg.output.xfm_ras,
    output:
        warp=temp(_lmk_xfm("warp.nii.gz")),
        invwarp=temp(_lmk_xfm("warp.nii.gz", from_="{template}", to="subject")),
        warped=temp(_lmk_warped("{reg}deformwarped")),
    log:
        _lmk_log("lmk_deform_reg", space="{template}", desc="{reg}"),
    group:
        "lmk_reg"
    params:
        iters=lambda wildcards: reg_iters(reg_presets[wildcards.reg], "deform_iters"),
        metric=lambda wildcards: reg_presets[wildcards.reg]["metric"],
        sigma1=lambda wildcards: reg_presets[wildcards.reg]["sigma1"],
        sigma2=lambda wildcards: reg_presets[wildcards.reg]["sigma2"],
        extra_pairs="",


use rule reg_ants as lmk_reg_ants with:
    input:
        fixed=get_template_for_reg,
        moving=_lmk_moving(),
        extra_fixed=[],
        extra_moving=[],
    output:
        xfm_itk=temp(_lmk_xfm("xfm.txt", desc="{reg}affine", type_="itk")),
        warp=temp(_lmk_xfm("warp.nii.gz")),
        invwarp=temp(_lmk_xfm("warp.nii.gz", from_="{template}", to="subject")),
        warped_affine=temp(_lmk_warped("{reg}affinewarped")),
        warped=temp(_lmk_warped("{reg}deformwarped")),
    log:
        _lmk_log("lmk_reg_ants", space="{template}", desc="{reg}"),
    params:
        args=lambda wildcards, input: ants_reg_args(
            reg_presets[wildcards.reg], input.fixed, input.moving, center_init=True
        ),
        identity=reg_identity_itk,


use rule reg_fireants as lmk_reg_fireants with:
    input:
        fixed=get_template_for_reg,
        moving=_lmk_moving(),
        extra_fixed=[],
        extra_moving=[],
    output:
        xfm_itk=temp(_lmk_xfm("xfm.txt", desc="{reg}affine", type_="itk")),
        warp=temp(_lmk_xfm("warp.nii.gz")),
        invwarp=temp(_lmk_xfm("warp.nii.gz", from_="{template}", to="subject")),
        warped_affine=temp(_lmk_warped("{reg}affinewarped")),
        warped=temp(_lmk_warped("{reg}deformwarped")),
    params:
        preset=lambda wildcards: {
            **reg_presets[wildcards.reg],
            "affine_iters": reg_iters(reg_presets[wildcards.reg], "affine_iters"),
            "deform_iters": reg_iters(reg_presets[wildcards.reg], "deform_iters"),
        },
        init_centers=lambda wildcards: reg_presets[wildcards.reg]["affine"],
        extra_fixed=[],
        extra_moving=[],


use rule compose_subject_to_template_warp as lmk_compose_subject_to_template_warp with:
    input:
        ref=rules.import_template_anat.output.anat,
        subject=_lmk_moving(),
        xfm_itk=_lmk_xfm("xfm.txt", desc="{reg}affine", type_="itk"),
        warp=_lmk_xfm("warp.nii.gz"),
        invwarp=_lmk_xfm("warp.nii.gz", from_="{template}", to="subject"),
    output:
        xfm_composite=lmk_xfm_fwd,
        xfm_composite_inv=lmk_xfm_inv,
    group:
        "lmk_reg"


# ---------------------------------------------------------------- MRI <-> SPIM

if lmk_do_mri:
    _mri_reg_wc = dict(
        dof="{dof}",
        iters="{iters}",
        radius="{radius}",
        gradsigma="{gradsigma}",
        warpsigma="{warpsigma}",
    )
    _mri_reg_cfg = {k: config["reg_mri"]["greedy"][k] for k in _mri_reg_wc}

    use rule affine_nlin_reg_mri_to_spim as lmk_affine_nlin_reg_mri_to_spim with:
        input:
            mri=bids(
                root=root,
                datatype="anat",
                desc="brain",
                suffix=f"{mri_suffix}.nii.gz",
                **inputs.subj_wildcards,
            ),
            spim=_lmk_moving(),
        output:
            xfm_ras=temp(
                _lmk_xfm(
                    "xfm.txt",
                    from_=mri_suffix,
                    to="SPIM",
                    desc="rigid",
                    type_="ras",
                    **_mri_reg_wc,
                )
            ),
            warp=temp(
                _lmk_xfm(
                    "warp.nii.gz",
                    from_=mri_suffix,
                    to="SPIM",
                    desc=None,
                    **_mri_reg_wc,
                )
            ),
            invwarp=temp(
                _lmk_xfm(
                    "warp.nii.gz",
                    from_="SPIM",
                    to=mri_suffix,
                    desc=None,
                    **_mri_reg_wc,
                )
            ),
            linwarped=temp(
                _lmk_warped(
                    "linearwarped",
                    space="SPIM",
                    suffix=f"{mri_suffix}.nii.gz",
                    **_mri_reg_wc,
                )
            ),
            warped=temp(
                _lmk_warped(
                    "deformwarped",
                    space="SPIM",
                    suffix=f"{mri_suffix}.nii.gz",
                    **_mri_reg_wc,
                )
            ),
        params:
            iters="{iters}",
            metric_rigid="NMI",
            dof="{dof}",
            metric="NCC {radius}",
            sigma1="{gradsigma}vox",
            sigma2="{warpsigma}vox",

    use rule compose_mri_to_spim_warp as lmk_compose_mri_to_spim_warp with:
        input:
            ref=_lmk_moving(),
            subject=bids(
                root=root,
                datatype="anat",
                desc="brain",
                suffix=f"{mri_suffix}.nii.gz",
                **inputs.subj_wildcards,
            ),
            xfm_ras=_lmk_xfm(
                "xfm.txt",
                from_=mri_suffix,
                to="SPIM",
                desc="rigid",
                type_="ras",
                **_mri_reg_cfg,
            ),
            warp=_lmk_xfm(
                "warp.nii.gz", from_=mri_suffix, to="SPIM", desc=None, **_mri_reg_cfg
            ),
            invwarp=_lmk_xfm(
                "warp.nii.gz", from_="SPIM", to=mri_suffix, desc=None, **_mri_reg_cfg
            ),
        output:
            xfm_composite=_lmk_xfm("xfm.nii.gz", from_=mri_suffix, to="SPIM", desc=None),
            xfm_composite_inv=_lmk_xfm(
                "xfm.nii.gz", from_="SPIM", to=mri_suffix, desc=None
            ),


# ---------------------------------------------------------------- TRE


def _lmk_tre_inputs(wildcards):
    wc = dict(wildcards.items())
    d = {
        "spim_landmarks": lmk_spim_file(wildcards),
        "template_landmarks": config["landmarks_template"],
    }
    mri = lmk_mri_file(wildcards)
    if mri:
        d["mri_landmarks"] = mri
    for s in lmk_stains_of(wildcards):
        for m in lmk_methods:
            f = dict(wc, regstain=s, reg=m, template=config["template"])
            d[f"fwd__{s}__{m}"] = _fwd.format(**f)
            d[f"inv__{s}__{m}"] = lmk_xfm_inv.format(**f)
        if mri:
            f = dict(wc, regstain=s)
            d[f"mri2spim__{s}"] = _lmk_xfm(
                "xfm.nii.gz", from_="SPIM", to=mri_suffix, desc=None
            ).format(**f)
            d[f"spim2mri__{s}"] = _lmk_xfm(
                "xfm.nii.gz", from_=mri_suffix, to="SPIM", desc=None
            ).format(**f)
    return d


rule landmark_tre:
    """TRE of the manual landmarks for every method, stain and path of a scan."""
    input:
        unpack(_lmk_tre_inputs),
    output:
        tsv=bids(
            root=lmk_root, datatype="tre", suffix="tre.tsv", **inputs["spim"].wildcards
        ),
    threads: 1
    resources:
        mem_mb=24000,
        runtime=120,
    params:
        scriptdir=os.path.join(workflow.basedir, "scripts"),
        stains=lmk_stains_of,
        methods=lmk_methods,
        mri_suffix=mri_suffix,
        template=config["template"],
    script:
        "../scripts/landmark_tre_subject.py"


def _lmk_scans():
    """SPIM scans (wildcard dicts) that have a landmarks file."""
    from types import SimpleNamespace

    zl = inputs["spim"].zip_lists
    keys = list(zl)
    out = []
    for vals in zip(*[zl[k] for k in keys]):
        w = dict(zip(keys, vals))
        if _lmk_find(
            config["landmarks_spim"], SimpleNamespace(**w), "micr", require="from-SPIM"
        ):
            out.append(w)
    return out


lmk_scans = _lmk_scans()
if not lmk_scans:
    raise ValueError(
        f"--landmark_analysis: no SPIM scan has a landmarks file in {config['landmarks_spim']}"
    )
_lmk_shared = LF.shared_files(
    [
        (
            str(w),
            _lmk_find(
                config["landmarks_spim"],
                SimpleNamespace(**w),
                "micr",
                require="from-SPIM",
            ),
        )
        for w in lmk_scans
    ]
)
if _lmk_shared:
    raise ValueError(
        "--landmark_analysis: a SPIM landmarks file matches several scans; its points are on one "
        "scan's grid, so name the scan in the file (ses-/sample-/acq-): "
        + "; ".join(f"{os.path.basename(f)} -> {s}" for f, s in _lmk_shared.items())
    )
print(f"Landmark analysis: {len(lmk_scans)} scans with SPIM landmarks")


def lmk_name_check():
    """Landmarks found on only one side of a comparison, per scan. Landmarks are
    matched by name, so a misspelt or missing one would otherwise be left out
    without notice. Rows: scan, landmark, found_in, missing_from."""
    import sys as _sys
    from types import SimpleNamespace

    _sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
    import landmark_tre as _LT

    tpl = set(_LT.load_landmarks(config["landmarks_template"]))
    rows = []
    for w in lmk_scans:
        ns = SimpleNamespace(**w)
        scan = "_".join(
            str(w[k]) for k in ("subject", "session", "sample", "acq") if k in w
        )
        spim = set(_LT.load_landmarks(lmk_spim_file(ns)))
        sets = [("SPIM", spim, "template", tpl)]
        mri_f = lmk_mri_file(ns) if lmk_do_mri else None
        if mri_f:
            sets.append(("MRI", set(_LT.load_landmarks(mri_f)), "SPIM", spim))
        elif lmk_do_mri:
            rows.append(
                dict(
                    scan=scan,
                    landmark="(all)",
                    found_in="",
                    missing_from="MRI (no MRI landmarks file)",
                )
            )
        for a, A, b, B in sets:
            rows += [
                dict(scan=scan, landmark=n, found_in=a, missing_from=b)
                for n in sorted(A - B)
            ]
            rows += [
                dict(scan=scan, landmark=n, found_in=b, missing_from=a)
                for n in sorted(B - A)
            ]
    return rows


lmk_unmatched = lmk_name_check()
if lmk_unmatched:
    print(
        f"Landmark analysis: WARNING {len(lmk_unmatched)} landmark(s) found on only one side are left out "
        "of the TRE (names must match); listed in landmarks/group/landmark_name_check.tsv and the report:"
    )
    for _r in lmk_unmatched[:20]:
        print(
            f"  {_r['scan']}: {_r['landmark']} in {_r['found_in'] or '-'}, not in {_r['missing_from']}"
        )
    if len(lmk_unmatched) > 20:
        print(f"  ... and {len(lmk_unmatched) - 20} more")
else:
    print(
        "Landmark analysis: every landmark name matches across SPIM, template"
        + (" and MRI" if lmk_do_mri else "")
    )

if lmk_do_mri:
    # use only the MRI scans the landmarks were placed on: the session (and any
    # acq/run/rec) named in each subject's MRI landmark file
    for _subj in sorted({w["subject"] for w in lmk_scans}):
        _f = lmk_mri_file(SimpleNamespace(subject=_subj))
        if not _f:
            continue
        _ents = _lmk_mri_entities(_f)
        if not _ents:
            continue
        try:
            _n = len(inputs["mri"].filter(subject=_subj, **_ents).expand())
        except Exception:
            _n = 0
        if _n == 0:
            raise ValueError(
                f"sub-{_subj}: no MRI input matches {_ents} from its landmark file "
                f"{os.path.basename(_f)}; check --filter-mri / --path-mri"
            )
        mri_subject_filters[_subj] = _ents
        print(
            f"Landmark analysis: sub-{_subj}: MRI from the landmark file: {_ents} ({_n} scan(s))"
        )

    # several MRI scans of a subject are averaged on the grid of the first one
    # (preproc_mri.smk), so the MRI landmarks have to be in that scan's space
    for _subj in sorted({w["subject"] for w in lmk_scans}):
        try:
            _mris = sorted(get_subject_mri(_subj).expand())
        except Exception:
            _mris = []
        if len(_mris) > 1:
            print(
                f"Landmark analysis: sub-{_subj}: {len(_mris)} MRI scans are averaged on the grid "
                f"of {os.path.basename(_mris[0])}; the MRI landmarks must be in that scan's space"
            )

lmk_tre_tsvs = [
    bids(
        root=lmk_root, datatype="tre", suffix="tre.tsv", **inputs["spim"].wildcards
    ).format(**w)
    for w in lmk_scans
]


rule landmark_markups:
    """Every warped landmark set of a scan as a 3D Slicer markups file (one per path x
    method x stain, the TRE of each point in its description), plus the target sets."""
    input:
        tsv=bids(
            root=lmk_root, datatype="tre", suffix="tre.tsv", **inputs["spim"].wildcards
        ),
    output:
        index=bids(
            root=lmk_root,
            datatype="markups",
            suffix="markups.tsv",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=10,
    params:
        template=config["template"],
        mri_suffix=mri_suffix,
    script:
        "../scripts/landmark_markups.py"


lmk_markups = [
    bids(
        root=lmk_root,
        datatype="markups",
        suffix="markups.tsv",
        **inputs["spim"].wildcards,
    ).format(**w)
    for w in lmk_scans
]


rule landmark_name_check:
    """Landmarks found on only one side of a comparison (left out of the TRE)."""
    input:
        spim=[lmk_spim_file(SimpleNamespace(**w)) for w in lmk_scans],
        mri=[f for f in (lmk_mri_file(SimpleNamespace(**w)) for w in lmk_scans) if f],
    output:
        tsv=os.path.join(lmk_root, "group", "landmark_name_check.tsv"),
    threads: 1
    resources:
        mem_mb=1000,
        runtime=5,
    params:
        template_landmarks=config["landmarks_template"],
    run:
        import pandas as pd

        cols = ["scan", "landmark", "found_in", "missing_from"]
        pd.DataFrame(lmk_name_check(), columns=cols).to_csv(
            output.tsv, sep="\t", index=False
        )


rule landmark_tre_summary:
    """All scans' TRE in one table, and median/IQR per path x method x stain (x landmark)."""
    input:
        tsvs=lmk_tre_tsvs,
    output:
        all_tsv=os.path.join(lmk_root, "group", "landmarks_tre_allscans.tsv"),
        summary=os.path.join(lmk_root, "group", "landmarks_tre_summary.tsv"),
        per_landmark=os.path.join(lmk_root, "group", "landmarks_tre_per_landmark.tsv"),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=15,
    script:
        "../scripts/landmark_tre_summary.py"


rule landmark_report:
    """Report: best method x stain, per-landmark spread, landmark-vs-registration
    classification, cross-subject consistency, L/R symmetry, MRI vs SPIM, flags with
    screenshots. Flags only; landmarks are never changed."""
    input:
        tsv=rules.landmark_tre_summary.output.all_tsv,
        name_check=rules.landmark_name_check.output.tsv,
        template_image=expand(
            rules.import_template_anat.output.anat, template=config["template"]
        ),
    output:
        report=os.path.join(lmk_root, "group", "report", "report.md"),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
    params:
        template_landmarks=config["landmarks_template"],
        run_root=root if root else ".",
        regions=config.get("landmark_regions"),
        template=config["template"],
    script:
        "../scripts/landmark_report.py"


rule all_landmarks:
    input:
        lmk_tre_tsvs,
        lmk_markups,
        rules.landmark_tre_summary.output,
        rules.landmark_report.output,
