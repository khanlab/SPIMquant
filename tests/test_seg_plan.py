"""Per-stain segmentation/object method resolution (seg_plan.py)."""

import sys
from pathlib import Path

import pytest


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


SCRIPTS_DIR = _find_repo_root(Path(__file__).parent) / "spimquant/workflow/scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from seg_plan import SegPlanError, parse_specs, resolve_seg_plan  # noqa: E402

SEG_DEFAULTS = {"default": "gmm+n3k1", "stains": {"Abeta": "lantern", "abeta": "lantern"}}
OBJ_DEFAULTS = {"default": "watershed", "stains": {"Abeta": "components", "abeta": "components"}}
STAINS = ["AutoF", "Abeta", "Iba1", "CD31"]
FOR_SEG = ["Abeta", "Iba1", "NeuN"]


def resolve(seg=None, obj=None, stains=STAINS, for_seg=FOR_SEG):
    plan, _ = resolve_seg_plan(stains, for_seg, seg, obj, SEG_DEFAULTS, OBJ_DEFAULTS)
    return plan


def test_defaults():
    plan = resolve()
    assert list(plan.stains) == ["Abeta", "Iba1"]
    assert plan["Abeta"].seg_methods == ["lantern"]
    assert plan["Abeta"].object_methods == ["components"]
    assert plan["Iba1"].seg_methods == ["gmm+n3k1"]
    assert plan["Iba1"].object_methods == ["watershed"]
    assert plan["Iba1"].seg_source == "default"


def test_bare_method_is_global_including_abeta():
    plan = resolve(seg=["otsu+k3i2"])
    assert plan["Abeta"].seg_methods == ["otsu+k3i2"]
    assert plan["Iba1"].seg_methods == ["otsu+k3i2"]
    assert plan["Abeta"].seg_source == "cli"


def test_stain_spec_beats_bare():
    plan = resolve(seg=["otsu+k3i2", "Abeta=lantern"])
    assert plan["Abeta"].seg_methods == ["lantern"]
    assert plan["Iba1"].seg_methods == ["otsu+k3i2"]


def test_comparisons_primary_first():
    plan = resolve(seg=["Abeta=lantern,otsu+k3i2"], obj=["Iba1=watershed,components"])
    assert plan["Abeta"].seg_method == "lantern"
    assert plan["Abeta"].pairs() == [("lantern", "components"), ("otsu+k3i2", "components")]
    assert plan["Iba1"].pairs() == [("gmm+n3k1", "watershed"), ("gmm+n3k1", "components")]
    assert plan["Iba1"].pairs(objects=False) == [("gmm+n3k1", None)]


def test_object_none():
    plan = resolve(obj=["Iba1=none"])
    assert not plan["Iba1"].has_objects
    assert plan["Iba1"].pairs() == []
    assert plan.stains_with_objects == ["Abeta"]


def test_named_stain_is_segmented():
    plan = resolve(seg=["CD31=th500"])
    assert "CD31" in plan
    plan = resolve(obj=["AutoF=none"])
    assert "AutoF" in plan


def test_named_stain_missing_from_data_is_a_note():
    plan, notes = resolve_seg_plan(
        STAINS, FOR_SEG, ["NeuN=th10"], None, SEG_DEFAULTS, OBJ_DEFAULTS
    )
    assert "NeuN" not in plan
    assert any("NeuN" in n for n in notes)


@pytest.mark.parametrize(
    "seg,obj,match",
    [
        (["gmm+n3k1+inst"], None, r"\+inst"),
        (["Abeta=lantern+inst"], None, r"\+inst"),
        (["foo"], None, "unknown method"),
        (None, ["cc"], "components"),
        (None, ["inst"], "watershed"),
        (None, ["Iba1=none,watershed"], "none"),
        (["Abeta=lantern", "Abeta=th5"], None, "more than once"),
        (["th5", "th6"], None, "more than one"),
        (["=th5"], None, "missing stain"),
    ],
)
def test_errors(seg, obj, match):
    with pytest.raises(SegPlanError, match=match):
        resolve(seg=seg, obj=obj)


def test_parse_specs_dedup():
    bare, by_stain = parse_specs(["Abeta=lantern,lantern,th5", "gmm+n2k2p5"], "--x")
    assert bare == ["gmm+n2k2p5"]
    assert by_stain == {"Abeta": ["lantern", "th5"]}


def test_seg_method_patterns():
    plan = resolve(seg=["Abeta=dog+s2p5r4k6", "Iba1=gmm+n2k2p5"])
    assert plan["Abeta"].seg_method == "dog+s2p5r4k6"


def test_provenance_marks_primary():
    plan = resolve(seg=["Abeta=lantern,otsu+k3i2"])
    prov = [p for p in plan.provenance() if p["Stain"] == "Abeta"]
    assert [p["Primary"] for p in prov] == [True, False]


# --- object distances (object_distances.py) -------------------------------

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from object_distances import add_object_distances, equivalent_radius  # noqa: E402


def test_object_distances():
    # two plaques on the x axis; one cell 1 mm from the first plaque
    df = pd.DataFrame(
        {
            "pos_x": [0.0, 10.0, 1.0],
            "pos_y": [0.0, 0.0, 0.0],
            "pos_z": [0.0, 0.0, 0.0],
            "volume": [4 / 3 * np.pi * 0.5**3, 4 / 3 * np.pi * 0.5**3, 0.001],
            "stain": ["Abeta", "Abeta", "Iba1"],
        }
    )
    out = add_object_distances(df, ["pos_x", "pos_y", "pos_z"])
    np.testing.assert_allclose(out["Abeta_dist"], [10.0, 10.0, 1.0])
    np.testing.assert_allclose(out["Abeta_edgedist"], [9.5, 9.5, 0.5])
    # the only Iba1 object has no other Iba1 object
    assert np.isnan(out.loc[2, "Iba1_dist"])
    np.testing.assert_allclose(out.loc[:1, "Iba1_dist"], [1.0, 9.0])
    np.testing.assert_allclose(equivalent_radius(4 / 3 * np.pi * 8), 2.0)


def test_object_distances_edge_clipped():
    df = pd.DataFrame(
        {
            "pos_x": [0.0, 0.2],
            "pos_y": [0.0, 0.0],
            "pos_z": [0.0, 0.0],
            "volume": [4 / 3 * np.pi, 0.0],
            "stain": ["Abeta", "Iba1"],
        }
    )
    out = add_object_distances(df, ["pos_x", "pos_y", "pos_z"])
    assert out.loc[1, "Abeta_edgedist"] == 0.0
