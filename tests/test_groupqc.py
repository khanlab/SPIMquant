"""Group QC report (groupqc_subject.py, groupqc_report.py) on a synthetic group."""

import json
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PIL")
sys.path.insert(
    0, str(Path(__file__).parents[1] / "spimquant" / "workflow" / "scripts")
)
from groupqc_report import build_data, metrics_tsv, render  # noqa: E402
from groupqc_subject import build_subject  # noqa: E402

SHAPE = (40, 48, 32)
AFFINE = np.diag([0.1, 0.1, 0.1, 1.0])
N_SLICES = {"coronal": 3, "axial": 2, "sagittal": 1}


def _save(path, data):
    nib.save(nib.Nifti1Image(np.asarray(data), AFFINE), str(path))
    return str(path)


def _brain():
    ijk = np.indices(SHAPE).astype(float)
    centre = (np.array(SHAPE) - 1)[:, None, None, None] / 2
    radii = np.array([16, 20, 12])[:, None, None, None]
    return (((ijk - centre) / radii) ** 2).sum(0) <= 1


@pytest.fixture
def group(tmp_path):
    rng = np.random.default_rng(0)
    brain = _brain()
    labels = np.where(brain, 1 + (np.indices(SHAPE)[1] > SHAPE[1] // 2), 0)
    tpl = {
        "mask": _save(tmp_path / "mask.nii.gz", brain.astype(np.uint8)),
        "anat": _save(
            tmp_path / "anat.nii.gz",
            np.where(brain, 100 + 50 * np.indices(SHAPE)[2], 0).astype(np.float32),
        ),
        "dseg": _save(tmp_path / "dseg.nii.gz", labels.astype(np.int16)),
    }
    jsons = []
    for i in range(7):
        sub = tmp_path / f"sub-{i:02d}"
        (sub / "qc").mkdir(parents=True)
        image = np.where(brain, 1000 + 500 * np.indices(SHAPE)[2], 0).astype(float)
        image += rng.normal(0, 50, SHAPE) * brain
        if i == 6:
            # the bad subject: half the brain is missing
            image[:, : SHAPE[1] // 2, :] = 0
        ff = np.where(brain, rng.uniform(0, 5, SHAPE), 0)
        table = pd.DataFrame(
            {
                "index": [0, 1, 2],
                "name": ["Unknown", "A", "B"],
                "Abeta+count": [5, 100 + i, 120],
                "Abeta+volume": [np.nan, 2.0, 3.0],
                "Abeta+fieldfrac": [0, 2.0, 3.0],
            }
        )
        table.to_csv(sub / "segstats.tsv", sep="\t", index=False)
        (sub / "qc" / f"sub-{i:02d}_stain-Abeta_desc-gmm_segslices.png").write_bytes(
            b""
        )
        result = build_subject(
            label=f"sub-{i:02d}",
            entities={"subject": f"{i:02d}"},
            template_mask=tpl["mask"],
            template_anat=tpl["anat"],
            stains=["AutoF", "Abeta"],
            warped=[
                _save(sub / "autof.nii.gz", image.astype(np.float32)),
                _save(sub / "abeta.nii.gz", (image / 2).astype(np.float32)),
            ],
            seg_layers=[{"stain": "Abeta", "desc": "gmm"}],
            fieldfrac=[_save(sub / "ff.nii.gz", ff.astype(np.float32))],
            n_slices=N_SLICES,
            segstats=str(sub / "segstats.tsv"),
            figure_prefix=str(sub / "qc" / f"sub-{i:02d}"),
        )
        path = sub / "groupqc.json"
        path.write_text(json.dumps(result))
        jsons.append(str(path))
    pd.DataFrame(
        {"participant_id": [f"sub-{i:02d}" for i in range(7)], "group": list("aabbaab")}
    ).to_csv(tmp_path / "participants.tsv", sep="\t", index=False)
    return tmp_path, tpl, jsons


def test_subject_json(group):
    _, _, jsons = group
    result = json.loads(Path(jsons[0]).read_text())
    assert set(result["layers"]) == {
        "stain:AutoF",
        "stain:Abeta",
        "fieldfrac:Abeta:gmm",
    }
    assert {p: len(v) for p, v in result["slices"].items()} == N_SLICES
    m = result["metrics"]
    assert m["reg_coverage"] > 99
    assert m["Abeta:count"] == 100 + 120  # the unknown label is left out
    assert m["Abeta:density"] == pytest.approx(220 / 5)
    assert m["Abeta:fieldfrac"] == pytest.approx((2 * 2 + 3 * 3) / 5)
    assert [f["label"] for f in result["figures"]] == ["Segmentation slices"]


def test_report_flags_bad_subject(group):
    tmp_path, tpl, jsons = group
    html_path = tmp_path / "group" / "qc" / "groupqc.html"
    html_path.parent.mkdir(parents=True)
    data = build_data(
        jsons,
        tpl["anat"],
        tpl["mask"],
        tpl["dseg"],
        str(html_path),
        title="synthetic",
        template_name="TPL",
        participants_tsv=str(tmp_path / "participants.tsv"),
        group_by=["group"],
    )
    by_label = {s["label"]: s for s in data["subjects"]}
    assert "reg_coverage" in by_label["sub-06"]["flags"]
    assert "reg_ncc_group" in by_label["sub-06"]["flags"]
    assert not by_label["sub-00"]["flags"]
    assert by_label["sub-00"]["cov"] == {"group": "a"}
    assert data["group_by"] == ["group"]
    # figures are linked relative to the report
    assert by_label["sub-00"]["figures"][0]["src"].startswith("../../sub-00/qc/")
    assert data["layers"][0]["key"] == "stain:AutoF"
    assert data["shared_scale"]["fieldfrac:Abeta:gmm"] > 0

    render(data, str(html_path))
    metrics_tsv(data, str(tmp_path / "groupqc.tsv"))
    page = html_path.read_text()
    assert '<script type="application/json" id="qc-data">' in page
    table = pd.read_csv(tmp_path / "groupqc.tsv", sep="\t")
    assert table.loc[table.subject == "sub-06", "n_flags"].item() >= 2
