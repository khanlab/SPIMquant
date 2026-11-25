"""Warped landmarks of one scan as 3D Slicer markups (landmark_markups rule).

From the scan's TRE table, writes next to the output index (a TSV listing them):

  <scan>_from-<A>_to-<B>_stain-<regstain>_desc-<method>_landmarks.mrk.json
      landmarks of space A carried to space B by that registration; each point's
      description holds its TRE
  <scan>_space-<B>_desc-target_landmarks.mrk.json
      the landmarks they are compared with, in space B

Spaces: SPIM (the level-5 lightsheet grid: open with the *_desc-N4brain_SPIM.nii.gz
of the scan), the template (tpl-<template> anat), and the MRI (desc-brain or
desc-preproc MRI of the subject). SPIM->template and MRI->template use the
fixed-point inverse, as the report does. Coordinates are LPS mm, as the inputs.

  python landmark_markups.py --tsv <scan>_tre.tsv --index <out>_markups.tsv \
      --template ABAv3 --mri-suffix T2w
"""
import argparse
import json
import os

import pandas as pd

SCHEMA = (
    "https://raw.githubusercontent.com/slicer/slicer/master/Modules/Loadable/"
    "Markups/Resources/Schema/markups-schema-v1.0.3.json#"
)
WARPED_RGB = [0.92, 0.41, 0.20]  # orange
TARGET_RGB = [0.11, 0.69, 0.48]  # aqua


def spaces(path, template, mri):
    """(from, to) spaces of a TRE path."""
    return {
        "spim_to_template": ("SPIM", template),
        "template_to_spim": (template, "SPIM"),
        "mri_to_spim": (mri, "SPIM"),
        "spim_to_mri": ("SPIM", mri),
        "mri_to_template": (mri, template),
    }[path]


def markups(points, name, rgb):
    """Slicer markups dict; points = [(label, description, (x, y, z))] in LPS mm."""
    return {
        "@schema": SCHEMA,
        "markups": [
            {
                "type": "Fiducial",
                "name": name,
                "coordinateSystem": "LPS",
                "coordinateUnits": "mm",
                "locked": True,
                "labelFormat": "%N-%d",
                "controlPoints": [
                    {
                        "id": str(i + 1),
                        "label": label,
                        "description": desc,
                        "position": [float(v) for v in xyz],
                        "orientation": [-1.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 1.0],
                        "selected": True,
                        "locked": True,
                        "visibility": True,
                        "positionStatus": "defined",
                    }
                    for i, (label, desc, xyz) in enumerate(points)
                ],
                "display": {
                    "visibility": True,
                    "opacity": 1.0,
                    "color": rgb,
                    "selectedColor": rgb,
                    "propertiesLabelVisibility": False,
                    "pointLabelsVisibility": True,
                    "textScale": 2.0,
                    "glyphType": "Sphere3D",
                    "glyphScale": 1.5,
                    "useGlyphScale": True,
                },
            }
        ],
    }


def write(obj, fn):
    with open(fn, "w") as fh:
        json.dump(obj, fh, indent=2)


def main(tsv, index, template, mri):
    df = pd.read_csv(tsv, sep="\t")
    df = df[df["inversion"].isin(["fixedpoint", "none"])]
    outdir = os.path.dirname(index) or "."
    os.makedirs(outdir, exist_ok=True)
    prefix = os.path.basename(index)[: -len("_markups.tsv")]
    rows, targets = [], {}
    for (path, method, stain), g in df.groupby(["path", "method", "stain"], sort=True):
        a, b = spaces(path, template, mri)
        pts = [
            (
                r.landmark,
                f"TRE {r.tre_um:.0f} um" + ("" if r.inside_field else " (outside the transform grid)"),
                (r.warped_x, r.warped_y, r.warped_z),
            )
            for r in g.itertuples()
        ]
        fn = f"{prefix}_from-{a}_to-{b}_stain-{stain}_desc-{method}_landmarks.mrk.json"
        write(markups(pts, f"{a}->{b} {method} {stain}", WARPED_RGB), os.path.join(outdir, fn))
        rows.append(dict(file=fn, path=path, from_space=a, to_space=b, method=method, stain=stain,
                         n=len(pts), median_tre_um=round(g["tre_um"].median(), 1)))
        for r in g.itertuples():
            targets.setdefault(b, {})[r.landmark] = (r.dst_x, r.dst_y, r.dst_z)
    for b, pts in sorted(targets.items()):
        fn = f"{prefix}_space-{b}_desc-target_landmarks.mrk.json"
        write(markups([(n, "target", xyz) for n, xyz in sorted(pts.items())], f"{b} target", TARGET_RGB),
              os.path.join(outdir, fn))
        rows.append(dict(file=fn, path="target", from_space="", to_space=b, method="", stain="",
                         n=len(pts), median_tre_um=""))
    pd.DataFrame(rows).to_csv(index, sep="\t", index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        main(sm.input.tsv, sm.output.index, sm.params.template, sm.params.mri_suffix)
    else:
        ap = argparse.ArgumentParser()
        ap.add_argument("--tsv", required=True)
        ap.add_argument("--index", required=True)
        ap.add_argument("--template", default="ABAv3")
        ap.add_argument("--mri-suffix", default="T2w")
        a = ap.parse_args()
        main(a.tsv, a.index, a.template, a.mri_suffix)
