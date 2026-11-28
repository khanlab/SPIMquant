"""Shared helpers of the group QC report (groupqc_subject.py, groupqc_report.py).

Every image of the report is a slice of a volume in template space, so the slices
of all subjects line up: the same slice positions, chosen from the template brain
mask, are used for every subject, which is what makes the montage and the flip
view comparable from one subject to the next.

Two kinds of layers are stored:

- gray: an intensity image, windowed per subject (robust percentiles and a gamma
  that puts the median at mid-gray, as in the registration QC report), so that
  differences in brightness between scans do not hide the anatomy.
- value: a quantitative map (field fraction), stored as 8-bit gray from 0 to the
  layer's ``scale`` and colored in the browser, so that the report can show it on
  each subject's own scale or on a scale shared by all subjects.
"""

import base64
from io import BytesIO

import numpy as np
from PIL import Image

# plane name, array axis through the plane (RAS+ arrays), direction of the axis
PLANES = [("coronal", 1, "A"), ("axial", 2, "S"), ("sagittal", 0, "R")]

# metric key -> label, section, unit, problem direction (low|high|both), help
METRICS = {
    "reg_ncc_template": (
        "NCC to template",
        "registration",
        "",
        "low",
        "Correlation of the warped registration stain with the template image, "
        "within the template brain mask.",
    ),
    "reg_ncc_group": (
        "NCC to group",
        "registration",
        "",
        "low",
        "Correlation of the warped registration stain with the mean of the other "
        "subjects, on the report slices: low when a subject is placed or shaped "
        "differently from the rest of the group.",
    ),
    "reg_coverage": (
        "Brain coverage",
        "registration",
        "%",
        "low",
        "Share of the template brain mask covered by image data: low when the "
        "tissue is cut, missing or outside the field of view.",
    ),
    "reg_jac_p1": (
        "Jacobian p1",
        "registration",
        "",
        "low",
        "1st percentile of the Jacobian determinant of the warp in the brain: "
        "strong local compression.",
    ),
    "reg_jac_p99": (
        "Jacobian p99",
        "registration",
        "",
        "high",
        "99th percentile of the Jacobian determinant of the warp in the brain: "
        "strong local expansion.",
    ),
    "reg_jac_folded": (
        "Folded voxels",
        "registration",
        "%",
        "high",
        "Share of brain voxels where the warp folds (Jacobian determinant at or "
        "below 0).",
    ),
}
# per-stain metrics: key suffix -> label, unit, direction, help
STAIN_METRICS = {
    "median": (
        "Median intensity",
        "",
        "both",
        "Median intensity in the template brain mask (raw units).",
    ),
    "p99": (
        "Intensity p99",
        "",
        "both",
        "99th percentile of intensity in the template brain mask (raw units).",
    ),
}
# per-segmentation metrics
SEG_METRICS = {
    "fieldfrac_mean": (
        "Mean field fraction",
        "%",
        "both",
        "Mean of the template-space field fraction map in the brain mask.",
    ),
    "count": (
        "Objects",
        "",
        "both",
        "Number of objects in all atlas regions (from the merged segstats table).",
    ),
    "density": (
        "Density",
        "/mm³",
        "both",
        "Objects per mm³ of all atlas regions.",
    ),
    "fieldfrac": (
        "Field fraction",
        "%",
        "both",
        "Volume-weighted field fraction over all atlas regions (from the merged "
        "segstats table).",
    ),
}


def slice_indices(mask, n_per_plane):
    """Slice indices for each plane, evenly spread through the mask."""
    out = {}
    for name, axis, _ in PLANES:
        n = n_per_plane.get(name, 0)
        if n <= 0:
            continue
        other = tuple(i for i in range(3) if i != axis)
        present = np.flatnonzero(mask.any(axis=other))
        if present.size == 0:
            present = np.arange(mask.shape[axis])
        out[name] = (
            np.linspace(present[0], present[-1], n + 2)[1:-1].round().astype(int)
        ).tolist()
    return out


def plane_axis(name):
    return next(axis for plane, axis, _ in PLANES if plane == name)


def to_view(data, axis, index):
    """Slice of a RAS+ array, with superior or anterior up (as in regqc_lib)."""
    return np.flipud(np.take(data, index, axis=axis).T)


def display_window(data, mask):
    """Window and gamma from the nonzero voxels in mask (as in regqc_lib)."""
    values = data[(data > 0) & mask]
    if values.size == 0:
        values = data[data > 0]
    if values.size == 0:
        return 0.0, 1.0, 1.0
    low, median, high = np.percentile(values, [1, 50, 99.5])
    high = max(high, low + 1e-6)
    scaled_median = (median - low) / (high - low)
    gamma = np.log(0.4) / np.log(scaled_median) if 0 < scaled_median < 1 else 1.0
    return float(low), float(high), float(np.clip(gamma, 0.2, 1.0))


def gray_u8(plane, window):
    low, high, gamma = window
    return (np.clip((plane - low) / (high - low), 0, 1) ** gamma * 255).astype(np.uint8)


def value_u8(plane, scale):
    return (np.clip(plane / scale, 0, 1) * 255).round().astype(np.uint8)


def uri(array_or_image, quality=80, lossless=False):
    """WebP data URI of a 2D uint8 array (gray) or a PIL image."""
    image = (
        array_or_image
        if isinstance(array_or_image, Image.Image)
        else Image.fromarray(array_or_image)
    )
    buffer = BytesIO()
    image.save(buffer, format="WEBP", quality=quality, lossless=lossless, method=6)
    return "data:image/webp;base64," + base64.b64encode(buffer.getvalue()).decode()


def decode_uri(data_uri):
    """2D uint8 array of a gray data URI."""
    raw = base64.b64decode(data_uri.split(",", 1)[1])
    return np.asarray(Image.open(BytesIO(raw)).convert("L"))


def outline_uri(labels, scale=2):
    """Transparent WebP with the boundaries between labels."""
    labels = labels.repeat(scale, axis=0).repeat(scale, axis=1)
    edges = np.zeros(labels.shape, dtype=bool)
    edges[:-1, :] |= labels[:-1, :] != labels[1:, :]
    edges[:, :-1] |= labels[:, :-1] != labels[:, 1:]
    rgba = np.zeros(labels.shape + (4,), dtype=np.uint8)
    rgba[edges] = (255, 214, 10, 235)
    return uri(Image.fromarray(rgba, "RGBA"), lossless=True)
