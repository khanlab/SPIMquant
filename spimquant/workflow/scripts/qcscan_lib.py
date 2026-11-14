"""Metrics and images of one scan for the dataset QC report (qc_scan.py).

Everything here works on small in-memory arrays: the downsampled volume of a
stain (in the axis order of the image), a tissue mask on the same grid, and 2D
full-resolution crops.

Images are stored as 8-bit with a square-root transfer of [0, hi], where hi is
recorded with the image, so that the report can decode them back to intensities
and apply the same display window to every scan of a stain.
"""

import base64
import io

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

PLANE_NAMES = {"L": "sagittal", "R": "sagittal", "A": "coronal", "P": "coronal"}
OUTLINE_RGB = (255, 214, 10)


def tissue_mask(volume):
    """Mask of the tissue in a downsampled volume.

    Otsu threshold of the smoothed log intensities, keeping the largest
    connected component with its holes filled.
    """
    from skimage.filters import threshold_otsu

    smoothed = ndi.gaussian_filter(np.log1p(np.clip(volume, 0, None)), 1.0)
    if smoothed.min() == smoothed.max():
        return np.zeros(volume.shape, bool)
    mask = smoothed > threshold_otsu(smoothed)
    labels, n = ndi.label(mask)
    if n > 1:
        sizes = ndi.sum_labels(mask, labels, index=np.arange(1, n + 1))
        mask = labels == (1 + int(np.argmax(sizes)))
    return ndi.binary_fill_holes(mask)


def edge_touch_fraction(mask):
    """Fraction of the tissue's surface voxels on the faces of the volume.

    High when the sample is cut off by the field of view.
    """
    if not mask.any():
        return 0.0
    surface = mask & ~ndi.binary_erosion(mask, border_value=0)
    on_face = np.zeros(mask.shape, bool)
    for axis in range(mask.ndim):
        index = [slice(None)] * mask.ndim
        for end in (0, -1):
            index[axis] = end
            on_face[tuple(index)] = True
    return float((surface & on_face).sum() / surface.sum())


def intensity_metrics(volume, mask):
    """Tissue and background intensities of a stain, and their contrast.

    snr is the tissue median over the background noise, and dynamic_range the
    brightest tissue (99.9th percentile) over the tissue median: high for sparse
    bright labelling on a dim background, near 1 for flat signal.
    """
    tissue = volume[mask]
    background = volume[~mask]
    if tissue.size == 0:
        return {}
    p1, p50, p99, p999 = (float(v) for v in np.percentile(tissue, [1, 50, 99, 99.9]))
    out = {"p1": p1, "p50": p50, "p99": p99, "p99_9": p999, "max": float(tissue.max())}
    out["dynamic_range"] = p999 / p50 if p50 > 0 else None
    if background.size:
        bg = float(np.median(background))
        noise = 1.4826 * float(np.median(np.abs(background - bg)))
        out["background"] = bg
        out["background_noise"] = noise
        out["tissue_to_background"] = p50 / bg if bg > 0 else None
        out["snr"] = (p50 - bg) / noise if noise > 0 else None
    return out


def axis_profiles(volume, mask, min_fraction=0.2):
    """Median tissue intensity along each axis, relative to the tissue median.

    Slices with little tissue (under min_fraction of the fullest slice) are
    None. Also returns the coefficient of variation of each profile: high when
    intensity falls off through the sample (attenuation, uneven clearing or
    illumination).
    """
    median = float(np.median(volume[mask])) if mask.any() else 0.0
    profiles, cvs = [], []
    for axis in range(volume.ndim):
        other = tuple(a for a in range(volume.ndim) if a != axis)
        counts = mask.sum(axis=other)
        keep = (
            counts >= max(1, min_fraction * counts.max())
            if counts.max()
            else counts > 0
        )
        masked = np.where(mask, volume, np.nan)
        moved = np.moveaxis(masked, axis, 0).reshape(volume.shape[axis], -1)
        profile = np.full(volume.shape[axis], np.nan)
        if keep.any() and median > 0:
            profile[keep] = np.nanmedian(moved[keep], axis=1) / median
        valid = profile[np.isfinite(profile)]
        cvs.append(float(valid.std() / valid.mean()) if valid.size > 2 else None)
        profiles.append([None if np.isnan(v) else round(float(v), 3) for v in profile])
    return profiles, cvs


def histogram(values, bins=64):
    """Counts in log-spaced bins of the positive values, for the report."""
    values = values[values > 0]
    if values.size == 0:
        return None
    lo = max(float(np.percentile(values, 0.1)), 1e-3)
    hi = float(values.max())
    if hi <= lo:
        return None
    counts, _ = np.histogram(values, bins=np.geomspace(lo, hi, bins + 1))
    return {"lo": lo, "hi": hi, "counts": [int(c) for c in counts]}


def correlation_matrix(volumes, mask):
    """Pearson correlation of the log intensities in tissue, between stains.

    Near 1 between two stains suggests bleed-through (or a shared background,
    such as autofluorescence).
    """
    logs = [np.log1p(np.clip(v[mask], 0, None)) for v in volumes]
    if not logs or logs[0].size < 3:
        return None
    corr = np.atleast_2d(np.corrcoef(np.stack(logs)))
    return [[None if np.isnan(v) else round(float(v), 3) for v in row] for row in corr]


def sharpness(image):
    """Fine detail relative to coarse structure in a crop, without units.

    The spread of the 1-2 pixel band over that of the 2-8 pixel band, so it does
    not depend on brightness; lower in an out-of-focus crop than in a sharp one
    of the same region. Compare it between scans of the same stain and region.
    """
    image = image.astype(np.float32)
    g1, g2, g8 = (ndi.gaussian_filter(image, s) for s in (1, 2, 8))
    coarse = float((g2 - g8).std())
    return float((g1 - g2).std() / coarse) if coarse > 0 else None


def saturated_fraction(image, ceiling):
    """Fraction of pixels at the top of the data type, or None if not known."""
    if ceiling is None or image.size == 0:
        return None
    return float((image >= ceiling).mean())


def slice_positions(mask, axis, n):
    """n slice indices along axis, evenly spaced through the tissue."""
    other = tuple(a for a in range(mask.ndim) if a != axis)
    present = np.flatnonzero(mask.any(axis=other))
    if present.size == 0:
        present = np.arange(mask.shape[axis])
    fractions = (np.arange(n) + 1) / (n + 1)
    return [int(round(present[0] + f * (present[-1] - present[0]))) for f in fractions]


def plane_name(orientation_letter):
    """Anatomical name of the plane perpendicular to an axis."""
    return PLANE_NAMES.get(orientation_letter, "axial")


def roi_center(labels, label_id):
    """Voxel deepest inside a region: the maximum of its distance transform.

    Unlike the centroid, it is always inside the region. None if the region is
    not in labels.
    """
    region = labels == label_id
    if not region.any():
        return None
    bounds = ndi.find_objects(region.astype(np.uint8))[0]
    padded = np.pad(region[bounds], 1)
    depth = ndi.distance_transform_edt(padded)
    local = np.unravel_index(int(np.argmax(depth)), depth.shape)
    return tuple(int(b.start + i - 1) for b, i in zip(bounds, local))


def spread_centers(mask, n):
    """n voxels spread through a mask, away from its surface and each other.

    The first is the deepest voxel; each next one is the voxel furthest from
    those already chosen, among voxels at least half as deep as the deepest.
    Deterministic, for crops of scans that are not registered to an atlas.
    """
    if not mask.any():
        return []
    depth = ndi.distance_transform_edt(mask)
    candidates = np.argwhere(depth >= 0.5 * depth.max())
    chosen = [np.array(np.unravel_index(int(np.argmax(depth)), mask.shape))]
    nearest = np.linalg.norm(candidates - chosen[0], axis=1)
    while len(chosen) < n and nearest.max() > 0:
        chosen.append(candidates[int(np.argmax(nearest))])
        nearest = np.minimum(nearest, np.linalg.norm(candidates - chosen[-1], axis=1))
    return [tuple(int(v) for v in c) for c in chosen]


def largest_labels(labels, n):
    """Ids of the n (or all, if None) labels with the most voxels, largest first."""
    ids, counts = np.unique(labels[labels > 0], return_counts=True)
    order = np.argsort(-counts, kind="stable")[:n]
    return [int(ids[i]) for i in order]


def hemisphere_free_name(name):
    """Name of a region without its hemisphere, e.g. L_Isocortex -> Isocortex."""
    lowered = name.lower()
    for prefix in ("l_", "r_", "left ", "right ", "left_", "right_"):
        if lowered.startswith(prefix):
            return name[len(prefix) :]
    for suffix in ("_l", "_r", " left", " right", "_left", "_right"):
        if lowered.endswith(suffix):
            return name[: -len(suffix)]
    return name


def match_labels(label_df, requested):
    """Rows of a label table matching names, abbreviations or indices."""
    rows = []
    for item in requested:
        text = str(item).lower()
        hit = None
        for _, row in label_df.iterrows():
            keys = [str(row["index"])] + [
                str(row[c]).lower() for c in ("name", "abbreviation") if c in row
            ]
            if text in keys:
                hit = row
                break
        if hit is None:
            raise ValueError(f"Atlas region {item!r} was not found in the label table")
        rows.append(hit)
    return rows


def resample_labels(labels, labels_affine, shape, affine):
    """Nearest-neighbour labels on the grid given by shape and affine."""
    to_labels = np.linalg.inv(labels_affine) @ affine
    grid = np.indices(shape).reshape(3, -1)
    coords = to_labels[:3, :3] @ grid + to_labels[:3, 3:4]
    out = ndi.map_coordinates(labels, coords, order=0, mode="constant", cval=0)
    return out.reshape(shape)


def display_slice(volume, axis, index, spacing, size, order=1):
    """A slice resampled to square pixels, with its longest side `size` pixels."""
    plane = np.take(volume, index, axis=axis)
    steps = [s for a, s in enumerate(spacing) if a != axis]
    extent = [n * s for n, s in zip(plane.shape, steps)]
    factor = size / max(extent)
    zoom = [s * factor for s in steps]
    return ndi.zoom(plane, zoom, order=order, mode="nearest", grid_mode=True)


def encode_image(image, hi, quality=82):
    """8-bit WebP data URI of an image, with a square-root transfer of [0, hi]."""
    scaled = np.sqrt(np.clip(image / hi, 0, 1)) if hi > 0 else np.zeros(image.shape)
    return _data_uri(
        Image.fromarray(np.round(scaled * 255).astype(np.uint8)),
        "WEBP",
        quality=quality,
        method=4,
    )


def encode_outline(labels):
    """Transparent PNG data URI of the boundaries between labels, or None."""
    edges = np.zeros(labels.shape, bool)
    edges[:-1] |= labels[:-1] != labels[1:]
    edges[:, :-1] |= labels[:, :-1] != labels[:, 1:]
    if not edges.any():
        return None
    rgba = np.zeros(labels.shape + (4,), np.uint8)
    rgba[edges] = OUTLINE_RGB + (255,)
    return _data_uri(Image.fromarray(rgba), "PNG", optimize=True)


def _data_uri(image, fmt, **kwargs):
    buffer = io.BytesIO()
    image.save(buffer, fmt, **kwargs)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/{fmt.lower()};base64,{encoded}"


def ncc(a, b, mask):
    """Normalized cross-correlation of two images within a mask."""
    if mask.sum() < 3:
        return None
    a, b = a[mask].astype(np.float64), b[mask].astype(np.float64)
    if a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def dice(a, b):
    total = a.sum() + b.sum()
    return float(2 * (a & b).sum() / total) if total else None
