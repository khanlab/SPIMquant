"""Landmark analysis helpers: which landmark file goes with which scan, which MRI
a landmark file names, convergence of the fixed-point inverse, empty TRE tables."""

import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

sys.path.insert(
    0, str(Path(__file__).parents[1] / "spimquant" / "workflow" / "scripts")
)
import landmark_files as LF  # noqa: E402
import landmark_tre as LT  # noqa: E402

A = "sub-01/ses-A/micr/sub-01_ses-A_sample-brain_acq-4x_from-SPIM_landmarks.json"
B = "sub-01/ses-B/micr/sub-01_ses-B_sample-brain_acq-4x_from-SPIM_landmarks.json"
NOSES = "sub-01/micr/sub-01_sample-brain_acq-4x_from-SPIM_landmarks.json"


def scan(**kw):
    return LF.scan_entities(SimpleNamespace(subject="01", **kw))


def test_scan_entities_short_and_long_names():
    assert scan(session="A", sample="brain", acq="4x") == {
        "ses": "A",
        "sample": "brain",
        "acq": "4x",
    }
    assert LF.scan_entities({"subject": "01", "acquisition": "4x"}) == {"acq": "4x"}


def test_file_of_another_session_is_not_used():
    # the bug: with one file, a ses-B scan used the ses-A landmarks
    assert (
        LF.pick_landmark_file([A], scan(session="B", sample="brain", acq="4x")) is None
    )
    assert LF.pick_landmark_file([A], scan(session="A", sample="brain", acq="4x")) == A


def test_each_session_gets_its_own_file():
    assert (
        LF.pick_landmark_file([A, B], scan(session="B", sample="brain", acq="4x")) == B
    )
    assert (
        LF.pick_landmark_file([A, B], scan(session="A", sample="brain", acq="4x")) == A
    )


def test_other_acquisition_or_sample_is_not_used():
    assert (
        LF.pick_landmark_file([A], scan(session="A", sample="brain", acq="9x")) is None
    )
    assert (
        LF.pick_landmark_file([A], scan(session="A", sample="hemi", acq="4x")) is None
    )


def test_file_without_a_tag_applies_and_named_file_is_preferred():
    assert (
        LF.pick_landmark_file([NOSES], scan(session="A", sample="brain", acq="4x"))
        == NOSES
    )
    assert (
        LF.pick_landmark_file([NOSES, A], scan(session="A", sample="brain", acq="4x"))
        == A
    )


def test_ambiguous_files_raise():
    a2 = A.replace("from-SPIM", "desc-v2_from-SPIM")
    with pytest.raises(ValueError):
        LF.pick_landmark_file([A, a2], scan(session="A", sample="brain", acq="4x"))


def test_shared_files():
    assert LF.shared_files([("s1", NOSES), ("s2", NOSES), ("s3", A), ("s4", None)]) == {
        NOSES: ["s1", "s2"]
    }


@pytest.mark.parametrize(
    "have,expected",
    [
        (
            ["subject", "session", "acquisition", "run", "reconstruction"],
            {"session": "11m", "acquisition": "GRE", "reconstruction": "x"},
        ),
        (
            ["subject", "session", "acq", "run", "rec"],
            {"session": "11m", "acq": "GRE", "rec": "x"},
        ),
        (["subject", "acquisition"], {"acquisition": "GRE"}),
    ],
)
def test_mri_entities_long_or_short_wildcards(have, expected):
    f = "sub-01/ses-11m/anat/sub-01_ses-11m_acq-GRE_rec-x_from-T2starw_landmarks.json"
    assert LF.mri_entities(f, have) == expected


def test_inv_converged():
    df = pd.DataFrame(
        dict(
            inversion=["fixedpoint", "fixedpoint", "shipped", "none", "fixedpoint"],
            inv_residual_um=[0.001, 264.0, None, None, LF.INV_TOL_UM],
        )
    )
    assert LF.inv_converged(df).tolist() == [True, False, True, True, True]
    assert LF.inv_converged(df.drop(columns="inv_residual_um")).all()


def test_empty_tre_table_has_a_header(tmp_path):
    f = tmp_path / "tre.tsv"
    LT.write_tsv([], f, fieldnames=["subject", "path", "tre_um", "inside_field"])
    df = pd.read_csv(f, sep="\t")
    assert df.empty and list(df.columns) == [
        "subject",
        "path",
        "tre_um",
        "inside_field",
    ]


@pytest.mark.parametrize(
    "label,name",
    [
        ("Template_CC_genu-1", "CC_genu"),
        ("SPIM_LV_anteromedial-1", "LV_anteromedial"),
        ("Template_SPIM_striatum_L_dorsomedial", "striatum_L_dorsomedial"),
        ("RH_rostral_tip-1", "RH_rostral_tip"),
        ("layer-5", "layer-5"),
        ("Template_layer-5-1", "layer-5"),
        ("lobule-10", "lobule-10"),
    ],
)
def test_landmark_name_keeps_own_numbers(label, name):
    assert LT.landmark_name(label) == name


def _rows(scan, path, method, stain, tre, n=3, inversion="none"):
    return [
        dict(
            scan=scan,
            path=path,
            inversion=inversion,
            method=method,
            stain=stain,
            tre_um=tre,
        )
        for _ in range(n)
    ]


def test_failed_template_registration_dropped_from_every_path():
    df = pd.DataFrame(
        _rows("s1", "spim_to_template", "ants", "PI", 5000, inversion="fixedpoint")
        + _rows(
            "s1", "template_to_spim", "ants", "PI", 1500
        )  # under the cut on its own
        + _rows("s1", "mri_to_template", "ants", "PI", 900, inversion="fixedpoint")
        + _rows("s1", "spim_to_template", "greedy", "PI", 150, inversion="fixedpoint")
        + _rows("s1", "mri_to_spim", "greedy", "PI", 200)
    )
    keep, failed = LF.failed_registrations(df, 2000)
    kept = df[keep]
    assert not ((kept.method == "ants") & (kept.stain == "PI")).any()
    assert set(kept.path) == {"spim_to_template", "mri_to_spim"}
    assert failed[["registration", "scan", "method", "stain"]].values.tolist() == [
        ["template", "s1", "ants", "PI"]
    ]
    assert failed["n"].tolist() == [9]


def test_failed_mri_registration_dropped_from_mri_paths():
    df = pd.DataFrame(
        _rows("s1", "mri_to_spim", "greedy", "PI", 4000)
        + _rows("s1", "spim_to_mri", "greedy", "PI", 300)
        + _rows("s1", "mri_to_template", "ants", "PI", 600, inversion="fixedpoint")
        + _rows("s1", "mri_to_template", "ants", "AutoF", 600, inversion="fixedpoint")
        + _rows("s1", "spim_to_template", "ants", "PI", 150, inversion="fixedpoint")
    )
    keep, failed = LF.failed_registrations(df, 2000)
    kept = df[keep]
    assert set(map(tuple, kept[["path", "stain"]].drop_duplicates().values)) == {
        ("mri_to_template", "AutoF"),
        ("spim_to_template", "PI"),
    }
    assert failed["registration"].tolist() == ["MRI"]


def test_nothing_failed():
    df = pd.DataFrame(
        _rows("s1", "spim_to_template", "ants", "PI", 150, inversion="fixedpoint")
    )
    keep, failed = LF.failed_registrations(df, 2000)
    assert keep.all() and failed.empty
