#!/usr/bin/env python3
"""Copy example outputs of a full-resolution run into docs/examples/.

    pixi run -e dev python docs/gen/collect_examples.py OUTPUT_DIR \
        --dataset "one APP/PS1 mouse brain, 4x SmartSPIM" \
        --resolution "1.6 x 1.6 x 2.75 um" \
        --command "spimquant bids out participant ..."

Use a run on a scan at about 2 x 2 x 4 um or finer: the test dataset is too
coarse to show what the reports look like. The script

- copies the reports and QC figures listed in ITEMS (PNGs downsized),
- renders a featuremap over the template,
- writes Markdown excerpts of the main tables,
- writes docs/examples/manifest.yml, from which docs/gen/build.py makes the
  Example outputs page,
- prints a table of run time and memory per rule from OUTPUT_DIR/benchmarks/
  (if the run used --benchmark), for docs/usage/running.md.

A required item that matches no file is an error: if an output was renamed,
update ITEMS rather than dropping it.
"""

import argparse
import datetime
import glob
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

DOCS = Path(__file__).resolve().parents[1]
EXAMPLES = DOCS / "examples"
ASSETS = EXAMPLES / "assets"
MAX_WIDTH = 1600

# (section, title, glob relative to the output dir, text, required)
ITEMS = [
    (
        "Registration",
        "Registration QC report",
        "sub-*/xfm/*_desc-*_regqc.html",
        "The subject warped to the template and the template warped to the "
        "subject, with atlas contours. Open this first for every brain.",
        True,
    ),
    (
        "Segmentation",
        "Segmentation overview: {stain}",
        "sub-*/qc/*_stain-*_desc-*_segslices.png",
        "Slices through the whole brain with the mask overlaid.",
        True,
    ),
    (
        "Segmentation",
        "Regions at full resolution: {stain}",
        "sub-*/qc/*_stain-*_seg-*_desc-*_roimontage.png",
        "Zoomed crops of atlas regions at the segmentation level, with the mask.",
        True,
    ),
    (
        "Segmentation",
        "Threshold: {stain}",
        "sub-*/seg/*_stain-*_desc-*_thresholds.png",
        "The intensity histogram of the corrected image, the fitted mixture and the threshold.",
        True,
    ),
    (
        "Segmentation",
        "Intensity histogram: {stain}",
        "sub-*/qc/*_stain-*_histogram.png",
        "Raw intensity histogram at full resolution.",
        True,
    ),
    (
        "Segmentation",
        "Objects: {stain}",
        "sub-*/qc/*_stain-*_obj-*_objectstats.png",
        "Distributions of object volume and shape.",
        True,
    ),
    (
        "Segmentation",
        "Depth profile: {stain}",
        "sub-*/qc/*_stain-*_desc-*_zprofile.png",
        "Field fraction against slice position, to spot depth-dependent artefacts.",
        True,
    ),
    (
        "Vessels",
        "Vessel overview: {stain}",
        "sub-*/qc/*_stain-*_desc-*_vesselslices.png",
        "Slices with the vesselFM mask overlaid.",
        False,
    ),
    (
        "Regions",
        "Per-region summary",
        "sub-*/qc/*_seg-*_from-*_roisummary.png",
        "Per-region values of every stain and metric for this brain.",
        True,
    ),
    (
        "Dataset QC",
        "Dataset QC report",
        "qc/qcreport.html",
        "The `qc` analysis level: every scan of the dataset on one page.",
        False,
    ),
    (
        "Group",
        "Group QC report",
        "group/qc/groupqc.html",
        "Registration and segmentation QC of all subjects on one page.",
        False,
    ),
]


def stain_of(path):
    m = re.search(r"_stain-([^_]+)_", Path(path).name)
    return m.group(1) if m else ""


def copy_asset(src, dest_dir):
    dest = dest_dir / Path(src).name
    if src.endswith(".png"):
        from PIL import Image

        im = Image.open(src)
        if im.width > MAX_WIDTH:
            im = im.resize((MAX_WIDTH, round(im.height * MAX_WIDTH / im.width)))
        im.save(dest, optimize=True)
    else:
        shutil.copy(src, dest)
    return dest


def render_featuremap(out_dir, dest_dir):
    """Abeta (or the first stain's) field fraction over the template."""
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt
    import nibabel as nib
    import numpy as np

    maps = sorted(
        glob.glob(f"{out_dir}/sub-*/featuremap/*_seg-coarse_*+fieldfrac.nii.gz")
    ) or sorted(glob.glob(f"{out_dir}/sub-*/featuremap/*+fieldfrac.nii.gz"))
    anat = glob.glob(f"{out_dir}/tpl-*/tpl-*_anat.nii.gz")
    if not maps or not anat:
        return None
    fmap, tpl = nib.load(maps[0]).get_fdata(), nib.load(anat[0]).get_fdata()
    if fmap.shape != tpl.shape:
        return None
    zs = np.linspace(0.2, 0.8, 6) * tpl.shape[1]
    fig, axes = plt.subplots(1, 6, figsize=(18, 3.6), facecolor="black")
    vmax = np.percentile(fmap[fmap > 0], 98) if (fmap > 0).any() else 1
    for ax, y in zip(axes, zs.astype(int)):
        ax.imshow(np.rot90(tpl[:, y, :]), cmap="gray")
        m = np.ma.masked_equal(np.rot90(fmap[:, y, :]), 0)
        im = ax.imshow(m, cmap="magma", vmin=0, vmax=vmax, alpha=0.85)
        ax.axis("off")
    cb = fig.colorbar(im, ax=axes, fraction=0.015, pad=0.01)
    cb.set_label("field fraction (%)", color="white")
    cb.ax.tick_params(colors="white")
    dest = dest_dir / "featuremap_fieldfrac.png"
    fig.savefig(dest, dpi=90, bbox_inches="tight", facecolor="black")
    return maps[0], dest


def table_excerpts(out_dir, dest_dir):
    import pandas as pd

    items = []
    merged = sorted(
        glob.glob(f"{out_dir}/sub-*/tabular/*_seg-coarse_from-*_mergedsegstats.tsv")
    ) or sorted(glob.glob(f"{out_dir}/sub-*/tabular/*_mergedsegstats.tsv"))
    if merged:
        df = pd.read_csv(merged[0], sep="\t")
        if "index" in df.columns:
            df = df[df["index"] != 0]  # label 0: outside the atlas
        metric_cols = [
            c
            for m in ("fieldfrac", "density", "count")
            for c in df.columns
            if c.endswith("+" + m) and not c.startswith("coloc")
        ][:6]
        cols = [c for c in ("index", "name") if c in df.columns] + metric_cols
        top = df.sort_values(metric_cols[0], ascending=False) if metric_cols else df
        md = dest_dir / "mergedsegstats.md"
        md.write_text(
            f"{len(df)} regions × {len(df.columns)} columns; the first metric columns, "
            f"regions sorted by `{metric_cols[0] if metric_cols else ''}`:\n\n"
            + top[cols].head(8).to_markdown(index=False, floatfmt=".3g")
            + "\n"
        )
        items.append(
            ("Tables", "Per-region table (`mergedsegstats.tsv`)", merged[0], md)
        )
    props = sorted(glob.glob(f"{out_dir}/sub-*/tabular/*_regionprops.parquet"))
    if props:
        df = pd.read_parquet(props[0])
        md = dest_dir / "regionprops.md"
        counts = df["stain"].value_counts().to_dict() if "stain" in df.columns else {}
        md.write_text(
            f"{len(df):,} objects ({', '.join(f'{k}: {v:,}' for k, v in counts.items())}) × "
            f"{len(df.columns)} columns: `{'`, `'.join(df.columns)}`.\n\n"
            + df.head(5).to_markdown(index=False, floatfmt=".4g")
            + "\n"
        )
        items.append(("Tables", "All objects (`regionprops.parquet`)", props[0], md))
    return items


def benchmark_table(out_dir):
    """Rule | jobs | wall time | peak memory, from Snakemake --benchmark files."""
    import pandas as pd

    rows = []
    for f in glob.glob(f"{out_dir}/benchmarks/**/*", recursive=True):
        if not Path(f).is_file():
            continue
        try:
            b = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if "s" not in b.columns:
            continue
        rule = (
            Path(f).parent.name
            if Path(f).parent != Path(out_dir) / "benchmarks"
            else Path(f).stem
        )
        rows.append(
            {
                "rule": rule,
                "s": b["s"].max(),
                "max_rss": b.get("max_rss", pd.Series([0])).max(),
            }
        )
    if not rows:
        return None
    df = (
        pd.DataFrame(rows)
        .groupby("rule")
        .agg(jobs=("s", "size"), s=("s", "max"), max_rss=("max_rss", "max"))
    )
    df = df.sort_values("s", ascending=False)
    lines = ["| Rule | Jobs | Longest job | Peak memory |", "|---|---|---|---|"]
    for rule, r in df.iterrows():
        if r.s < 60:
            continue
        t = f"{r.s / 3600:.1f} h" if r.s >= 3600 else f"{r.s / 60:.0f} min"
        lines.append(f"| `{rule}` | {int(r.jobs)} | {t} | {r.max_rss / 1024:.1f} GB |")
    lines.append(
        f"\nTotal over all jobs: {df.s.mul(df.jobs).sum() / 3600:.1f} h of wall time (jobs run in parallel)."
    )
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("output_dir")
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--resolution", required=True)
    ap.add_argument("--stains", default="")
    ap.add_argument("--command", required=True)
    args = ap.parse_args()
    out = Path(args.output_dir)

    if ASSETS.exists():
        shutil.rmtree(ASSETS)
    ASSETS.mkdir(parents=True)
    entries, missing = [], []
    for section, title, pattern, text, required in ITEMS:
        files = sorted(glob.glob(str(out / pattern)))
        files = [f for f in files if "nomask" not in f]
        if not files:
            if required:
                missing.append(pattern)
            continue
        seen = set()
        for f in files:
            stain = stain_of(f)
            if "{stain}" in title and stain in seen:
                continue  # one per stain
            seen.add(stain)
            if not any(
                e["section"] == section and e["title"] == title.format(stain=stain)
                for e in entries
            ):
                dest = copy_asset(f, ASSETS)
                entries.append(
                    {
                        "section": section,
                        "title": title.format(stain=stain),
                        "text": text,
                        "file": f"assets/{dest.name}",
                        "source": str(Path(f).relative_to(out)),
                    }
                )
            if "{stain}" not in title:
                break
    if missing:
        sys.exit("no file matches required pattern(s):\n  " + "\n  ".join(missing))

    fm = render_featuremap(out, ASSETS)
    if fm:
        entries.append(
            {
                "section": "Regions",
                "title": "Field fraction in the template",
                "text": "Each region's field fraction painted into the template "
                "(`featuremap/`), over the template image.",
                "file": f"assets/{fm[1].name}",
                "source": str(Path(fm[0]).relative_to(out)),
            }
        )
    for section, title, src, md in table_excerpts(out, ASSETS):
        entries.append(
            {
                "section": section,
                "title": title,
                "file": f"assets/{md.name}",
                "source": str(Path(src).relative_to(out)),
            }
        )

    order = [
        "Registration",
        "Segmentation",
        "Vessels",
        "Regions",
        "Tables",
        "Dataset QC",
        "Group",
    ]
    entries.sort(key=lambda e: order.index(e["section"]))
    version = subprocess.run(
        ["git", "describe", "--tags", "--always"],
        cwd=DOCS.parent,
        capture_output=True,
        text=True,
    ).stdout.strip()
    manifest = {
        "intro": "Outputs of one full-resolution brain, as SPIMquant writes them. "
        "Every file is described in the [output reference](../reference/outputs.md).",
        "dataset": args.dataset,
        "resolution": args.resolution,
        "stains": args.stains,
        "spimquant_version": version,
        "command": args.command,
        "collected": datetime.date.today().isoformat(),
        "items": entries,
    }
    (EXAMPLES / "manifest.yml").write_text(
        yaml.safe_dump(manifest, sort_keys=False, width=100)
    )
    size = sum(p.stat().st_size for p in ASSETS.iterdir()) / 1e6
    print(f"wrote {len(entries)} items, {size:.1f} MB, to {EXAMPLES}")
    table = benchmark_table(out)
    if table:
        print(
            "\nResource table for docs/usage/running.md (between the RESOURCES-TABLE markers):\n"
        )
        print(table)


if __name__ == "__main__":
    main()
