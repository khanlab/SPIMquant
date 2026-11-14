"""QC metrics and images of one scan, for the dataset QC report.

Writes one JSON with, for every stain of the scan, metrics and overview slices
from the downsampled image and full-resolution crops, as embedded images (see
qcscan_lib.py for their encoding). The crops are at the same places in every
stain: in atlas regions if the scan was registered (rule qc_register), or else
spread through the tissue.
"""

import json

import nibabel as nib
import numpy as np
import pandas as pd
from dask_setup import get_dask_client
from zarrnii import ZarrNii

import qcscan_lib as lib


def _affine(znimg):
    return np.asarray(znimg.affine.matrix, dtype=float)


def _physical(affine, index):
    return tuple(float(v) for v in (affine @ np.append(index, 1.0))[:3])


def _crop_regions(inp, params, shape, affine, mask):
    """Places to crop: (name, abbreviation, voxel on the downsampled grid or None).

    Also returns the atlas labels on the downsampled grid, or None."""
    if "dseg" not in inp.keys():
        centers = lib.spread_centers(mask, params.n_crops)
        return [
            (f"Crop {i + 1}", f"crop{i + 1}", c) for i, c in enumerate(centers)
        ], None

    dseg_nii = nib.load(inp.dseg)
    labels = lib.resample_labels(
        np.asarray(dseg_nii.dataobj), dseg_nii.affine, shape, affine
    )
    label_df = pd.read_csv(inp.label_tsv, sep="\t")
    if params.crop_labels:
        rows = lib.match_labels(label_df, params.crop_labels)
    else:
        # the largest regions of the template, so that every scan gets the same
        # (one hemisphere of each, in atlases that split them)
        by_id = {int(row["index"]): row for _, row in label_df.iterrows()}
        template_labels = np.asarray(nib.load(inp.template_dseg).dataobj)
        rows, seen = [], set()
        for i in lib.largest_labels(template_labels, None):
            row = by_id.get(i, {"index": i})
            base = lib.hemisphere_free_name(str(row.get("name", i)))
            if base not in seen:
                seen.add(base)
                rows.append(row)
            if len(rows) == params.n_crops:
                break
    regions = []
    for row in rows:
        label_id = int(row["index"])
        name = str(row.get("name", label_id))
        abbreviation = str(row.get("abbreviation", name))
        regions.append((name, abbreviation, lib.roi_center(labels, label_id)))
    return regions, labels


def _registration(inp, params, mask, labels):
    """How well the quick registration did, or None if there was none."""
    if labels is None:
        return None
    # every 4th voxel is plenty for a correlation, and the template can be large
    step = (slice(None, None, 4),) * 3
    fixed = np.asarray(nib.load(inp.fixed).dataobj[step], dtype=np.float32)
    warped = np.asarray(nib.load(inp.warped).dataobj[step], dtype=np.float32)
    return {
        "stain": params.reg_stain,
        "template": params.template,
        "atlas": params.atlas_seg,
        # similarity to the template, and overlap of the atlas with the tissue
        "ncc": lib.ncc(fixed, warped, fixed > 0),
        "dice": lib.dice(labels > 0, mask),
    }


def main():
    inp, params = snakemake.input, snakemake.params
    stains = list(params.stains)
    quality = params.image_quality
    ceiling = (2 ** (8 * params.itemsize) - 1) if params.itemsize <= 2 else None

    low = ZarrNii.from_file(
        inp.spim,
        level=params.level,
        channel_labels=stains,
        downsample_near_isotropic=True,
        **params.zarrnii_kwargs,
    )
    volumes = np.asarray(low.data.compute(), dtype=np.float32)
    affine = _affine(low)
    spacing = [float(np.linalg.norm(affine[:3, a])) for a in range(3)]
    orientation = str(low.xyz_orientation)
    # data axes are in the image's axis order; orientation is given for x, y, z
    axis_letters = [
        orientation["xyz".index(d)] for d in low.dims if d in ("x", "y", "z")
    ]
    shape = volumes.shape[1:]

    # tissue from the registration stain, or else all stains together
    if params.reg_stain in stains:
        mask = lib.tissue_mask(volumes[stains.index(params.reg_stain)])
    else:
        scale = np.percentile(volumes, 99, axis=(1, 2, 3)).clip(1e-6)
        mask = lib.tissue_mask((volumes / scale[:, None, None, None]).mean(axis=0))

    regions, labels = _crop_regions(inp, params, shape, affine, mask)

    # overview slices, and atlas outlines on them
    planes = []
    for axis in range(3):
        positions = lib.slice_positions(mask, axis, params.overview_slices)
        outlines = [
            (
                lib.encode_outline(
                    lib.display_slice(
                        labels, axis, i, spacing, params.overview_size, order=0
                    )
                )
                if labels is not None
                else None
            )
            for i in positions
        ]
        planes.append(
            {
                "name": lib.plane_name(axis_letters[axis]),
                "axis": axis,
                "positions": positions,
                "outlines": outlines,
            }
        )

    # the same order of planes in every scan, whatever its orientation
    planes.sort(key=lambda p: p["name"])

    # full-resolution crops: one plane of the image at each place
    high = ZarrNii.from_file(
        inp.spim,
        level=params.crop_level,
        channel_labels=stains,
        **params.zarrnii_kwargs,
    )
    high_shape = dict(zip(high.dims, high.data.shape))
    patch = (
        min(params.crop_size, high_shape["x"]),
        min(params.crop_size, high_shape["y"]),
        1,
    )
    high_affine = _affine(high)
    crop_spacing = [float(np.linalg.norm(high_affine[:3, a])) for a in range(3)]
    crops = []
    for _, _, center in regions:
        if center is None:
            crops.append(None)
            continue
        cropped = high.crop_centered(_physical(affine, center), patch_size=patch)
        data = np.asarray(cropped.data.compute(), dtype=np.float32)
        crops.append(data.squeeze(axis=list(cropped.dims).index("z")))

    channels = {}
    for c, stain in enumerate(stains):
        volume = volumes[c]
        metrics = lib.intensity_metrics(volume, mask)
        profiles, cvs = lib.axis_profiles(volume, mask)
        valid_cvs = [v for v in cvs if v is not None]
        metrics["profile_cv"] = max(valid_cvs) if valid_cvs else None

        stain_crops = [None if crop is None else crop[c] for crop in crops]
        present = [crop for crop in stain_crops if crop is not None]
        # one transfer for all images of the stain in this scan
        peaks = [float(volume[mask].max()) if mask.any() else float(volume.max())]
        peaks += [float(np.percentile(crop, 99.99)) for crop in present]
        hi = max(max(peaks), 1e-6)

        crop_items = []
        for crop in stain_crops:
            if crop is None:
                crop_items.append(None)
                continue
            crop_items.append(
                {
                    "image": lib.encode_image(crop, hi, quality),
                    "sharpness": lib.sharpness(crop),
                    "saturated": lib.saturated_fraction(crop, ceiling),
                }
            )
        sharp = [i["sharpness"] for i in crop_items if i and i["sharpness"] is not None]
        sat = [i["saturated"] for i in crop_items if i and i["saturated"] is not None]
        metrics["sharpness"] = float(np.median(sharp)) if sharp else None
        metrics["saturated"] = float(np.max(sat)) if sat else None

        channels[stain] = {
            "metrics": metrics,
            "hi": hi,
            # default display window, in intensities
            "window": [metrics.get("p1", 0.0), metrics.get("p99_9", hi)],
            "histogram": lib.histogram(volume[mask]),
            "profiles": [profiles[plane["axis"]] for plane in planes],
            "overview": [
                [
                    lib.encode_image(
                        lib.display_slice(
                            volume, plane["axis"], i, spacing, params.overview_size
                        ),
                        hi,
                        quality,
                    )
                    for i in plane["positions"]
                ]
                for plane in planes
            ],
            "crops": crop_items,
        }

    voxel_mm3 = float(np.prod(spacing))
    entities = {k: v for k, v in params.entities.items() if v}
    result = {
        "entities": entities,
        "stains": stains,
        "registration": _registration(inp, params, mask, labels),
        "meta": {
            "orientation": orientation,
            "shape": [int(high_shape[d]) for d in ("x", "y", "z")],
            "voxel_um": [
                round(1000 * crop_spacing[list(high.dims[-3:]).index(d)], 3)
                for d in ("x", "y", "z")
            ],
            "crop_level": int(params.crop_level),
            "overview_level": int(params.level),
            "overview_voxel_um": [round(1000 * s, 2) for s in spacing],
            "bits": 8 * int(params.itemsize),
        },
        "scan": {
            "tissue_volume_mm3": float(mask.sum()) * voxel_mm3,
            "edge_touch": lib.edge_touch_fraction(mask),
            "correlation": lib.correlation_matrix(list(volumes), mask),
        },
        "planes": [{"name": p["name"], "outlines": p["outlines"]} for p in planes],
        "regions": [
            {"name": name, "abbreviation": abbreviation}
            for name, abbreviation, _ in regions
        ],
        "channels": channels,
    }
    with open(snakemake.output.json, "w") as f:
        json.dump(result, f, separators=(",", ":"))


if __name__ == "__main__":
    with get_dask_client("threads", snakemake.threads):
        main()
