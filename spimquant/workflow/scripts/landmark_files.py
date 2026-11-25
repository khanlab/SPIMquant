"""Helpers of the landmark analysis that need no images: which landmark file
belongs to which scan, which MRI a landmark file names, and which fixed-point
inversions converged. Used by rules/landmarks.smk and the landmark scripts, and
tested in tests/test_landmark_files.py."""

import os
import re

# a fixed-point inverse whose residual is above this did not converge; its rows
# are kept in the tables but left out of the summaries and the report
INV_TOL_UM = 1.0

# BIDS tag in a file name -> the names a snakebids input may use for it (the
# wildcard is the tag, acq/rec, but the entity name is acquisition/reconstruction)
_MRI_TAGS = {
    "ses": ("session", "ses"),
    "acq": ("acquisition", "acq"),
    "run": ("run",),
    "rec": ("reconstruction", "rec"),
}


def name_entities(path):
    """{tag: value} of the BIDS entities in a file name (sub, ses, acq, ...)."""
    out = {}
    for part in os.path.basename(path).split("_"):
        k, sep, v = part.partition("-")
        if sep and v and re.fullmatch(r"[a-zA-Z]+", k):
            out[k] = v
    return out


def scan_entities(wildcards):
    """{tag: value} of the session, sample and acquisition of a SPIM scan, from
    rule wildcards or a zip_lists row (either the short or the long names)."""
    get = wildcards.get if isinstance(wildcards, dict) else lambda k: getattr(wildcards, k, None)
    out = {}
    for tag, names in (("ses", ("session", "ses")), ("sample", ("sample",)), ("acq", ("acq", "acquisition"))):
        for n in names:
            v = get(n)
            if v not in (None, ""):
                out[tag] = str(v)
                break
    return out


def pick_landmark_file(paths, scan):
    """The landmark file of a scan among a subject's landmark files, or None.

    scan: {tag: value} of the scan (scan_entities). A file naming a different
    value for one of those tags (e.g. ses-B for a ses-A scan) belongs to another
    scan and is dropped. A file that does not name a tag applies to every scan;
    if several files remain, the ones naming the scan's value are preferred.
    Raises when it stays ambiguous."""
    hits = [p for p in sorted(paths) if all(name_entities(p).get(t, v) == v for t, v in scan.items())]
    for t, v in scan.items():
        if len(hits) > 1:
            hits = [p for p in hits if name_entities(p).get(t) == v] or hits
    if len(hits) > 1:
        raise ValueError(f"several landmark files for scan {scan}: {hits}")
    return hits[0] if hits else None


def shared_files(scan_files):
    """{file: [scan, ...]} for landmark files picked for more than one scan.
    scan_files: [(scan label, file)]. Points sit on one scan's grid, so a file
    used by two scans would give the other one a wrong TRE."""
    by_file = {}
    for scan, f in scan_files:
        if f:
            by_file.setdefault(f, []).append(scan)
    return {f: s for f, s in by_file.items() if len(s) > 1}


def mri_entities(path, have):
    """{MRI input entity: value} named in an MRI landmark file name (ses, acq,
    run, rec), keyed by whichever name the MRI input uses (have: its wildcard
    names, e.g. acq or acquisition). Tags the input lacks are left out."""
    have = set(have)
    out = {}
    for tag, v in name_entities(path).items():
        for n in _MRI_TAGS.get(tag, ()):
            if n in have:
                out[n] = v
                break
    return out


def inv_converged(df, tol_um=INV_TOL_UM):
    """Boolean Series: False for fixed-point rows whose inverse residual is above
    tol_um (the inversion did not converge), True otherwise."""
    import pandas as pd

    if "inv_residual_um" not in df.columns:
        return pd.Series(True, index=df.index)
    res = pd.to_numeric(df["inv_residual_um"], errors="coerce")
    return ~((df["inversion"] == "fixedpoint") & (res > tol_um))
