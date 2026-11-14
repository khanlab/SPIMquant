"""Dataset QC report (the qc analysis level): qcscan_lib.py and qcreport_lib.py."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PIL")
sys.path.insert(
    0, str(Path(__file__).parents[1] / "spimquant" / "workflow" / "scripts")
)
import qcreport_lib as report  # noqa: E402
import qcscan_lib as lib  # noqa: E402


def _phantom(shape=(40, 48, 56), seed=0):
    """A bright ball of tissue on a dim noisy background."""
    rng = np.random.default_rng(seed)
    grid = np.indices(shape) - np.array(shape)[:, None, None, None] / 2
    ball = np.sqrt((grid**2).sum(axis=0)) < 14
    volume = rng.normal(100, 5, shape) + 900 * ball
    return volume.astype(np.float32), ball


def test_tissue_mask_and_metrics():
    volume, ball = _phantom()
    mask = lib.tissue_mask(volume)
    assert lib.dice(mask, ball) > 0.9
    assert lib.edge_touch_fraction(mask) == 0
    metrics = lib.intensity_metrics(volume, mask)
    assert metrics["p50"] == pytest.approx(1000, rel=0.05)
    assert metrics["tissue_to_background"] == pytest.approx(10, rel=0.1)
    assert metrics["snr"] > 50


def test_edge_touch_when_cut_off():
    mask = np.zeros((20, 20, 20), bool)
    mask[:10, 5:15, 5:15] = True
    assert lib.edge_touch_fraction(mask) > 0.1


def test_profiles_show_attenuation():
    volume, ball = _phantom()
    ramp = np.linspace(1.5, 0.5, volume.shape[0])[:, None, None]
    _, flat = lib.axis_profiles(volume, ball)
    _, sloped = lib.axis_profiles(volume * ramp, ball)
    assert sloped[0] > 3 * flat[0]
    assert sloped[0] == max(sloped)


def test_sharpness_drops_with_blur():
    from scipy import ndimage as ndi

    rng = np.random.default_rng(1)
    image = ndi.gaussian_filter(rng.normal(size=(128, 128)), 1.0) + 10
    assert lib.sharpness(ndi.gaussian_filter(image, 3)) < 0.5 * lib.sharpness(image)
    # and does not depend on brightness
    assert lib.sharpness(5 * image) == pytest.approx(lib.sharpness(image), rel=1e-3)


def test_crop_centers():
    labels = np.zeros((20, 20, 20), int)
    labels[2:8, 2:8, 2:8] = 3
    # an L shape, whose centroid is outside it
    labels[10:18, 10:12, 5] = 4
    labels[10:12, 10:18, 5] = 4
    assert labels[lib.roi_center(labels, 3)] == 3
    assert labels[lib.roi_center(labels, 4)] == 4
    assert lib.roi_center(labels, 9) is None
    assert lib.largest_labels(labels, 1) == [3]
    assert lib.largest_labels(labels, None) == [3, 4]

    _, ball = _phantom()
    centers = lib.spread_centers(ball, 5)
    assert len(set(centers)) == 5
    assert all(ball[c] for c in centers)
    assert centers == lib.spread_centers(ball, 5)


def test_hemisphere_free_name():
    assert lib.hemisphere_free_name("L_Isocortex") == "Isocortex"
    assert lib.hemisphere_free_name("Right Thalamus") == "Thalamus"
    assert lib.hemisphere_free_name("Lateral septum") == "Lateral septum"


def test_resample_labels_between_axis_orders():
    labels = np.arange(2 * 3 * 4).reshape(2, 3, 4)
    # the same grid stored with reversed axes (x, y, z <-> z, y, x)
    flip = np.eye(4)[[2, 1, 0, 3]]
    out = lib.resample_labels(labels, np.eye(4), (4, 3, 2), flip)
    assert np.array_equal(out, labels.transpose(2, 1, 0))


def test_encode_image_roundtrip():
    import base64
    import io

    from PIL import Image

    image = np.linspace(0, 1000, 64 * 64, dtype=np.float32).reshape(64, 64)
    uri = lib.encode_image(image, 1000, quality=100)
    stored = np.asarray(Image.open(io.BytesIO(base64.b64decode(uri.split(",")[1]))))
    if stored.ndim == 3:
        stored = stored[..., 0]
    decoded = 1000 * (stored.astype(float) / 255) ** 2
    assert np.abs(decoded - image).max() < 40


def test_canonical_stains():
    mapping = report.canonical_stains(
        ["YoPro", "YOPRO", "abeta", "Abeta", "Iba1"], {"YOPRO": ["YoPro", "YoPro1"]}
    )
    assert mapping["YoPro"] == mapping["YOPRO"] == "YOPRO"
    assert mapping["abeta"] == mapping["Abeta"] == "abeta"
    assert mapping["Iba1"] == "Iba1"


def test_flags():
    others = [10, 11, 9, 10, 10.5, 2]
    assert report.flag(2, others, "low", "snr")
    assert report.flag(2, others, "high", "snr") is None
    assert report.flag(10, others, "both", "snr") is None
    # too few scans to call anything an outlier
    assert report.flag(2, [10, 11, 2], "low", "snr") is None
    # fixed limits apply whatever the number of scans
    assert report.flag(0.01, [0.01], "high", "saturated")
    assert report.flag(None, others, "low", "snr") is None


def _scan(subject, stains, snr, registered=True):
    image = lib.encode_image(np.ones((4, 4), np.float32), 1.0)
    return {
        "entities": {"subject": subject, "sample": "brain"},
        "stains": stains,
        "registration": (
            {"stain": stains[0], "template": "T", "atlas": "a", "ncc": 0.5, "dice": 0.9}
            if registered
            else None
        ),
        "meta": {
            "orientation": "RAS",
            "shape": [4, 4, 4],
            "voxel_um": [1, 1, 1],
            "crop_level": 0,
            "overview_level": 1,
            "overview_voxel_um": [2, 2, 2],
            "bits": 16,
        },
        "scan": {"tissue_volume_mm3": 400.0, "edge_touch": 0.0, "correlation": None},
        "planes": [
            {"name": n, "outlines": [None]} for n in ("axial", "coronal", "sagittal")
        ],
        "regions": [{"name": "Isocortex", "abbreviation": "CTX"}],
        "channels": {
            s: {
                "metrics": {"snr": snr, "p50": 100.0, "saturated": 0.0},
                "hi": 1.0,
                "window": [0.0, 1.0],
                "histogram": None,
                "profiles": [[1.0], [1.0], [1.0]],
                "overview": [[image]] * 3,
                "crops": [{"image": image, "sharpness": 0.5, "saturated": 0.0}],
            }
            for s in stains
        },
    }


def test_report_with_different_stains(tmp_path):
    scans = [
        _scan("01", ["YoPro", "Abeta"], 10.0),
        _scan("02", ["YOPRO", "Abeta"], 9.8),
        _scan("03", ["YoPro", "Iba1"], 10.1),
        _scan("04", ["YoPro", "Abeta"], 10.2),
        _scan("05", ["YoPro", "Abeta"], 0.5),
        _scan("06", ["NeuN"], 10.0, registered=False),
    ]
    participants = pd.DataFrame(
        {
            "participant_id": [f"sub-0{i}" for i in range(1, 7)],
            "protocol": list("aabbcc"),
        }
    )
    data = report.build_report_data(
        scans, "test", {"YOPRO": ["YoPro"]}, participants, None
    )
    assert data["stains"] == ["YOPRO", "Abeta", "Iba1", "NeuN"]
    by_id = {s["id"]: s for s in data["scans"]}
    assert set(by_id["sub-03_sample-brain"]["channels"]) == {"YOPRO", "Iba1"}
    assert by_id["sub-03_sample-brain"]["label"] == "sub-03"
    assert by_id["sub-03_sample-brain"]["groups"] == {"protocol": "b"}
    # flagged against the other scans of the same stain only
    assert "snr" in by_id["sub-05_sample-brain"]["channels"]["YOPRO"]["flags"]
    assert not by_id["sub-01_sample-brain"]["channels"]["YOPRO"]["flags"]
    # Iba1 and NeuN have one scan each: nothing to compare them to
    assert not by_id["sub-06_sample-brain"]["channels"]["NeuN"]["flags"]

    table = report.metrics_table(data)
    assert len(table) == sum(len(s["stains"]) for s in scans)
    assert set(table["stain"]) == set(data["stains"])

    out = tmp_path / "qcreport.html"
    report.render_report(out, data)
    page = out.read_text()
    assert "sub-06" in page and "NeuN" in page
    # self-contained: nothing loaded from elsewhere
    assert 'src="http' not in page and 'href="http' not in page
