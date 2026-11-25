"""Landmark report for --landmark_analysis: which tool/stain registers best, and
what the landmarks themselves look like. Flags points for a human look; never
changes them.

Run by the landmark_report rule, or:
  python landmark_report.py --tsv landmarks_tre_allscans.tsv --template-landmarks DIR \
      --template-image tpl-ABAv3_anat.nii.gz --run-root OUTPUT_DIR --out REPORT_DIR

Sections
 summary figure (methods, paths, subjects, regions)
 0. failed registrations (listed, left out of the rest)
 1. best method x stain per path (median TRE heatmaps)
 2. per-landmark TRE spread (boxplots)
 3. TRE by region (--landmark_regions, or a built-in map)
 4. TRE by subject
 5. landmark vs registration: per subject x landmark, the consensus of all
    method x stain warps vs the landmark. If independent registrations agree with
    each other but not with the landmark, the landmark is the likely problem.
 6. consistency across subjects in template space (cluster spread, bias, outliers)
 7. left/right symmetry of paired landmarks
 8. MRI vs lightsheet agreement
 9. flagged points with screenshots
Also writes markups/: the subjects' consensus points and the template points in
template space, as 3D Slicer markups.
"""
import argparse
import glob
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PRIMARY = ("spim_to_template", "fixedpoint")
# a registration whose median TRE on a scan exceeds this did not converge (e.g. a
# preset without an affine stage on unaligned input); reported, not used in the stats
FAIL_UM = 2000.0

# default landmark -> region map (first match wins); --landmark_regions overrides
DEFAULT_REGIONS = [
    ("cerebellum", r"^cerebellum"),
    ("lateral cortex", r"^(visual)?cortex_"),
    ("hemisphere surface", r"^[LR]H_"),
    ("ventricles / aqueduct", r"^([LR]V|3V|CA)_"),
    ("corpus callosum", r"^CC_"),
    ("striatum", r"^striatum"),
    ("thalamus", r"^thalamus"),
    ("midbrain", r"^midbrain"),
]
PATHS = [  # (path, inversion, label) in report order
    ("spim_to_template", "fixedpoint", "LS → template"),
    ("template_to_spim", "none", "template → LS"),
    ("mri_to_spim", "none", "MRI → LS"),
    ("spim_to_mri", "none", "LS → MRI"),
    ("mri_to_template", "fixedpoint", "MRI → template"),
]
# figure colours (categorical slots in a fixed order; text stays ink)
C_LS, C_MRI, C_L, C_R = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True, "font.size": 9,
    "axes.titlesize": 10, "axes.titleweight": "bold",
})


def load_regions(path):
    if not path:
        return DEFAULT_REGIONS
    t = pd.read_csv(path, sep="\t")
    return list(zip(t["region"].astype(str), t["pattern"].astype(str)))


def region_of(name, regions):
    for region, pat in regions:
        if re.search(pat, name):
            return region
    return "other"


def write_markups(points, fn, name, rgb):
    """[(label, description, (x, y, z) LPS mm)] -> Slicer markups file."""
    import json

    cps = [dict(id=str(i + 1), label=l, description=d, position=[float(v) for v in p],
                selected=True, locked=True, visibility=True, positionStatus="defined")
           for i, (l, d, p) in enumerate(points)]
    js = {"@schema": "https://raw.githubusercontent.com/slicer/slicer/master/Modules/Loadable/Markups/"
                     "Resources/Schema/markups-schema-v1.0.3.json#",
          "markups": [dict(type="Fiducial", name=name, coordinateSystem="LPS", coordinateUnits="mm",
                           locked=True, controlPoints=cps,
                           display=dict(color=rgb, selectedColor=rgb, pointLabelsVisibility=True,
                                        textScale=2.0, glyphType="Sphere3D", glyphScale=1.5))]}
    with open(fn, "w") as fh:
        json.dump(js, fh, indent=2)


def mad(x):
    x = np.asarray(x, float)
    return np.median(np.abs(x - np.median(x))) * 1.4826


def xyz(df, prefix):
    return df[[f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"]].to_numpy(float)


def pair_of(name):
    """(base, side) for paired landmark names, else None."""
    for pat, rep in [(r"_L_", "_{}_"), (r"_R_", "_{}_"), (r"_L$", "_{}"), (r"_R$", "_{}")]:
        if re.search(pat, name):
            side = "L" if "L" in pat else "R"
            return re.sub(pat, rep, name, count=1), side
    m = re.match(r"^(L|R)(H|V)_(.+)$", name)
    if m:
        return "{}" + m.group(2) + "_" + m.group(3), m.group(1)
    return None


def md_table(df, floatfmt="{:.0f}"):
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append(
            "| "
            + " | ".join(floatfmt.format(v) if isinstance(v, (float, np.floating)) else str(v) for v in r)
            + " |"
        )
    return "\n".join(out)


def scan_id(df):
    keys = [k for k in ("subject", "session", "sample", "acq") if k in df.columns]
    return df[keys].astype(str).agg("_".join, axis=1)


def load_landmarks(path):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import landmark_tre as LT

    return LT.load_landmarks(path)


# ------------------------------------------------------------------ screenshots
def _slice_png(img_path, point_lps, ax_row, title, marks=()):
    import nibabel as nib

    img = nib.load(img_path)
    data = np.asarray(img.dataobj, dtype=np.float32)
    if data.ndim > 3:
        data = data[..., 0]
    inv = np.linalg.inv(img.affine)
    to_vox = lambda p: (inv @ np.r_[np.asarray(p) * [-1, -1, 1], 1])[:3]
    v = to_vox(point_lps)
    c = np.clip(np.round(v).astype(int), 0, np.array(data.shape) - 1)
    half = 40
    patch = data[tuple(slice(max(ci - half, 0), ci + half) for ci in c)]
    pos = patch[patch > 0]
    vmin, vmax = (np.percentile(pos, [1, 99]) if pos.size else (0, 1))
    for k, ax in enumerate(ax_row):
        sl = [slice(None)] * 3
        sl[k] = c[k]
        im = data[tuple(sl)]
        rows = [d for d in range(3) if d != k]
        ax.imshow(im.T, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)
        ax.plot(v[rows[0]], v[rows[1]], "r+", ms=14, mew=2)
        for col, p in marks:
            u = to_vox(p)
            ax.plot(u[rows[0]], u[rows[1]], "o", mfc="none", mec=col, ms=10, mew=1.5)
        ax.set_xlim(c[rows[0]] - half, c[rows[0]] + half)
        ax.set_ylim(c[rows[1]] - half, c[rows[1]] + half)
        ax.set_xticks([])
        ax.set_yticks([])
    ax_row[0].set_ylabel(title, fontsize=8)


def main(tsv, tpl_dir, tpl_image, run_root, out, regions_file=None, template="template"):
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    for old in glob.glob(os.path.join(out, "figures", "flag_*.png")):  # from an earlier report
        os.remove(old)
    regions = load_regions(regions_file)
    df = pd.read_csv(tsv, sep="\t")
    df["scan"] = scan_id(df)
    df = df[df["inside_field"].astype(bool)].copy()
    reg_key = ["scan", "path", "inversion", "method", "stain"]
    med = df.groupby(reg_key)["tre_um"].transform("median")
    failed = (df[med > FAIL_UM].groupby(reg_key)["tre_um"].agg(["median", "count"])
              .rename(columns={"median": "median_um", "count": "n"}).reset_index())
    df = df[med <= FAIL_UM].copy()
    tpl = load_landmarks(tpl_dir)
    md = ["# Landmark registration report", ""]
    nsub = df["subject"].nunique()
    md += [
        f"Scans: {df['scan'].nunique()} ({nsub} subjects). Methods: {', '.join(sorted(df['method'].unique()))}. "
        f"Stains: {', '.join(sorted(df['stain'].unique()))}. Landmarks: {df['landmark'].nunique()}.",
        "",
        "TRE in µm. Primary measure: SPIM→template with the fixed-point inverse. Rows with a point "
        "outside a transform's grid are excluded. **Landmarks are ground truth here and are never "
        "edited; flags are for a human to look at.**",
        "",
    ]
    header_end = len(md)

    if not failed.empty:
        failed.to_csv(os.path.join(out, "failed_registrations.tsv"), sep="\t", index=False)
        fs = (failed.groupby(["method", "stain"])["median_um"].median().reset_index()
              .rename(columns={"median_um": "median_TRE_um"}))
        md += ["## 0. Failed registrations", "",
               f"Median TRE above {FAIL_UM:.0f} µm on a scan: the registration did not converge. "
               "These are left out of every section below (all rows are in "
               "`failed_registrations.tsv` and the group TSVs).", "",
               md_table(fs), ""]

    # 1. method x stain
    md += ["## 1. Which method and stain register best", ""]
    best_rows = []
    for (path, inv), g in df.groupby(["path", "inversion"]):
        piv = g.groupby(["method", "stain"])["tre_um"].median().unstack()
        if piv.size > 1:
            fig, ax = plt.subplots(figsize=(1.2 + 1.1 * piv.shape[1], 0.8 + 0.45 * piv.shape[0]))
            imh = ax.imshow(piv.values, cmap="viridis_r")
            ax.set_xticks(range(piv.shape[1]), piv.columns)
            ax.set_yticks(range(piv.shape[0]), piv.index)
            for i in range(piv.shape[0]):
                for j in range(piv.shape[1]):
                    ax.text(j, i, f"{piv.values[i, j]:.0f}", ha="center", va="center",
                            color="w" if piv.values[i, j] > np.nanmedian(piv.values) else "k", fontsize=8)
            ax.set_title(f"median TRE (µm): {path} [{inv}]", fontsize=9)
            fig.colorbar(imh, ax=ax, fraction=0.04)
            fig.tight_layout()
            fn = f"figures/heatmap_{path}_{inv}.png"
            fig.savefig(os.path.join(out, fn), dpi=110)
            plt.close(fig)
            md += [f"![{path} {inv}]({fn})", ""]
        s = g.groupby(["method", "stain"])["tre_um"].agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75), "count"])
        s.columns = ["median", "q25", "q75", "n"]
        s = s.sort_values("median").reset_index()
        b = s.iloc[0]
        best_rows.append(dict(path=path, inversion=inv, best_method=b["method"], best_stain=b["stain"],
                              median_um=b["median"], q25=b["q25"], q75=b["q75"],
                              worst=f"{s.iloc[-1]['method']}/{s.iloc[-1]['stain']} ({s.iloc[-1]['median']:.0f})"))
    md += [md_table(pd.DataFrame(best_rows)), ""]
    pd.DataFrame(best_rows).to_csv(os.path.join(out, "best_method_stain.tsv"), sep="\t", index=False)
    best = {(r["path"], r["inversion"]): (r["best_method"], r["best_stain"]) for r in best_rows}
    paths = [(p, i, lab) for p, i, lab in PATHS if (p, i) in best]

    def best_rows_of(path, inv):
        m, st = best[(path, inv)]
        return df[(df.path == path) & (df.inversion == inv) & (df.method == m) & (df.stain == st)]

    df["region"] = df["landmark"].map(lambda n: region_of(n, regions))

    # 2. per-landmark boxplots
    md += ["## 2. TRE per landmark", ""]
    for path, inv in [PRIMARY, ("mri_to_spim", "none"), ("mri_to_template", "fixedpoint")]:
        g = df[(df.path == path) & (df.inversion == inv)]
        if g.empty:
            continue
        order = g.groupby("landmark")["tre_um"].median().sort_values().index
        fig, ax = plt.subplots(figsize=(10, 0.28 * len(order) + 1.2))
        ax.boxplot([g[g.landmark == n]["tre_um"] for n in order], vert=False, tick_labels=order,
                   flierprops=dict(ms=2))
        ax.set_xlabel("TRE (µm), all subjects x methods x stains")
        ax.set_title(f"{path} [{inv}]", fontsize=9)
        ax.tick_params(axis="y", labelsize=7)
        fig.tight_layout()
        fn = f"figures/box_{path}_{inv}.png"
        fig.savefig(os.path.join(out, fn), dpi=110)
        plt.close(fig)
        md += [f"![{path}]({fn})", ""]

    # 3. regions
    md += ["## 3. TRE by region", ""]
    reg_rows = []
    for p, i, lab in paths:
        g = best_rows_of(p, i)
        for region, h in g.groupby("region"):
            reg_rows.append(dict(region=region, path=lab, median_um=h["tre_um"].median(),
                                 n_landmarks=h["landmark"].nunique()))
    region_order = []
    if reg_rows:
        rt = pd.DataFrame(reg_rows)
        rt.to_csv(os.path.join(out, "regions.tsv"), sep="\t", index=False)
        wide = rt.pivot(index="region", columns="path", values="median_um")
        wide = wide[[lab for _, _, lab in paths if lab in wide.columns]]
        nl = rt.groupby("region")["n_landmarks"].max()
        region_order = list(wide.sort_values(wide.columns[0]).index)
        wide = wide.loc[region_order]
        wide.insert(0, "landmarks", nl.loc[region_order].astype(int))
        g = best_rows_of(*PRIMARY)
        m, st = best[PRIMARY]
        fig, ax = plt.subplots(figsize=(7, 0.42 * len(region_order) + 1.2))
        data = [g[g.region == r]["tre_um"].to_numpy() for r in region_order]
        bp = ax.boxplot(data, vert=False, tick_labels=region_order, widths=0.55, patch_artist=True,
                        medianprops=dict(color=INK, lw=1.5), flierprops=dict(ms=2, mec=INK2),
                        whiskerprops=dict(color=INK2), capprops=dict(color=INK2))
        for b in bp["boxes"]:
            b.set(facecolor=C_LS + "55", edgecolor=C_LS)
        ax.set_xlabel("TRE (µm), every landmark of every subject")
        ax.set_title(f"LS → template by region ({m}, {st})", loc="left")
        ax.grid(axis="y", visible=False)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "figures/regions.png"), dpi=130)
        plt.close(fig)
        md += [
            "Median TRE per region with each path's best method and stain (section 1). Regions come "
            "from " + (f"`{os.path.basename(regions_file)}`" if regions_file else "the built-in map")
            + "; a landmark goes to the first region whose pattern matches its name.",
            "",
            md_table(wide.reset_index()),
            "",
            "![regions](figures/regions.png)",
            "",
        ]

    # 4. subjects
    md += ["## 4. TRE by subject", ""]
    subj_rows = []
    for p, i, lab in paths:
        for subj, h in best_rows_of(p, i).groupby("subject"):
            subj_rows.append(dict(subject=subj, path=lab, median_um=h["tre_um"].median()))
    subj_order = []
    if subj_rows:
        stt = pd.DataFrame(subj_rows)
        stt.to_csv(os.path.join(out, "subjects.tsv"), sep="\t", index=False)
        wide = stt.pivot(index="subject", columns="path", values="median_um")
        wide = wide[[lab for _, _, lab in paths if lab in wide.columns]]
        subj_order = list(wide.sort_values(wide.columns[0]).index)
        g = best_rows_of(*PRIMARY)
        gm = best_rows_of("mri_to_template", "fixedpoint") if ("mri_to_template", "fixedpoint") in best else None
        fig, ax = plt.subplots(figsize=(7, 0.42 * len(subj_order) + 1.4))
        y = np.arange(len(subj_order))
        for k, (gg, col, lab) in enumerate([(g, C_LS, "LS → template"), (gm, C_MRI, "MRI → template")]):
            if gg is None or gg.empty:
                continue
            off = (-0.17 if gm is not None and not gm.empty else 0) * (1 if k == 0 else -1)
            data = [gg[gg.subject == s_]["tre_um"].to_numpy() for s_ in subj_order]
            bp = ax.boxplot(data, positions=y + off, vert=False, widths=0.3, patch_artist=True,
                            medianprops=dict(color=INK, lw=1.5), flierprops=dict(ms=2, mec=INK2),
                            whiskerprops=dict(color=INK2), capprops=dict(color=INK2), manage_ticks=False)
            for b in bp["boxes"]:
                b.set(facecolor=col + "55", edgecolor=col)
            ax.plot([], [], "s", color=col, label=lab)
        ax.set_yticks(y, [f"sub-{s_}" for s_ in subj_order])
        ax.set_xlabel("TRE (µm), every landmark (best method and stain of each path)")
        ax.set_title("TRE by subject", loc="left")
        ax.grid(axis="y", visible=False)
        ax.legend(frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "figures/subjects.png"), dpi=130)
        plt.close(fig)
        md += [
            "Median TRE per subject with each path's best method and stain. A subject much worse than "
            "the others on every path points at that scan (image quality, landmark set); worse only "
            "on the MRI paths, at its MRI.",
            "",
            md_table(wide.loc[subj_order].reset_index()),
            "",
            "![subjects](figures/subjects.png)",
            "",
        ]

    # 5. landmark vs registration (consensus)
    g = df[(df.path == PRIMARY[0]) & (df.inversion == PRIMARY[1])]
    cons = []
    for (scan, subj, name), h in g.groupby(["scan", "subject", "landmark"]):
        W = xyz(h, "warped")
        c = np.median(W, axis=0)
        t = xyz(h, "dst")[0]
        cons.append(dict(scan=scan, subject=subj, landmark=name, n_regs=len(h),
                         consensus_err_um=np.linalg.norm(c - t) * 1000,
                         method_spread_um=np.median(np.linalg.norm(W - c, axis=1)) * 1000,
                         best_tre_um=h["tre_um"].min(),
                         cx=c[0], cy=c[1], cz=c[2], tx=t[0], ty=t[1], tz=t[2]))
    cons = pd.DataFrame(cons)
    cons["ratio"] = cons["consensus_err_um"] / cons["method_spread_um"].clip(lower=10)
    per = cons.groupby("landmark").agg(
        median_consensus_err_um=("consensus_err_um", "median"),
        median_method_spread_um=("method_spread_um", "median"),
        median_ratio=("ratio", "median"),
        subjects=("subject", "nunique"),
    ).reset_index()
    hi_err = per["median_consensus_err_um"].quantile(0.75)
    hi_spread = per["median_method_spread_um"].quantile(0.75)

    def klass(r):
        if r.median_consensus_err_um >= hi_err and r.median_ratio >= 2:
            return "consistently hard (registrations agree, landmark disagrees: placement or ambiguous anatomy)"
        if r.median_method_spread_um >= hi_spread:
            return "tool/stain-dependent (registrations disagree)"
        return "ok"

    per["class"] = per.apply(klass, axis=1)
    per = per.sort_values("median_consensus_err_um", ascending=False)
    per.to_csv(os.path.join(out, "landmark_classes.tsv"), sep="\t", index=False)
    md += [
        "## 5. Landmark or registration?",
        "",
        "For each subject x landmark: *consensus* = median warped position over all method x stain "
        "registrations; *consensus error* = its distance to the template landmark; *method spread* = "
        "median distance of the registrations to the consensus. Independent registrations rarely make "
        "the same error, so a large consensus error with a small spread points at the landmark "
        "(placement or ambiguous anatomy), a large spread at the registration.",
        "",
        md_table(per[["landmark", "median_consensus_err_um", "median_method_spread_um", "median_ratio", "class"]]),
        "",
    ]

    # 4. template-space consistency across subjects
    md += ["## 6. Consistency across subjects (template space)", ""]
    rows4, outl = [], []
    if cons["subject"].nunique() >= 3:
        for name, h in cons.groupby("landmark"):
            C = h[["cx", "cy", "cz"]].to_numpy()
            t = h[["tx", "ty", "tz"]].to_numpy()[0]
            m = C.mean(0)
            d = np.linalg.norm(C - m, axis=1) * 1000
            rows4.append(dict(landmark=name, n=len(h), spread_rms_um=np.sqrt((d ** 2).mean()),
                              bias_of_mean_um=np.linalg.norm(m - t) * 1000))
            # each subject against the mean of the others (leave-one-out), so an
            # outlier does not pull the reference towards itself
            n = len(C)
            loo = (C.sum(0) - C) / (n - 1)
            d_loo = np.linalg.norm(C - loo, axis=1) * 1000
            thr = max(np.median(d_loo) + 3 * mad(d_loo), 150)
            for (_, r), di in zip(h.iterrows(), d_loo):
                if di > thr:
                    outl.append(dict(subject=r.subject, landmark=name, dist_to_others_um=di, threshold_um=thr))
        t4 = pd.DataFrame(rows4).sort_values("spread_rms_um", ascending=False)
        md += [
            "Spread = RMS distance of the subjects' consensus points to their mean; bias = distance of "
            "that mean from the template landmark (a large bias with small spread suggests the template "
            "point itself is off).",
            "",
            md_table(t4),
            "",
        ]
        # scatter
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        for k, (a, b) in enumerate([(0, 1), (0, 2), (1, 2)]):
            ax = axs[k]
            for name, h in cons.groupby("landmark"):
                C = h[["cx", "cy", "cz"]].to_numpy()
                t = h[["tx", "ty", "tz"]].to_numpy()[0]
                ax.scatter(C[:, a], C[:, b], s=6)
                ax.plot(t[a], t[b], "k+", ms=7)
            ax.set_xlabel("LPS " + "xyz"[a] + " (mm)")
            ax.set_ylabel("LPS " + "xyz"[b] + " (mm)")
            ax.set_aspect("equal")
        fig.suptitle("Consensus warped landmarks per subject (dots) vs template landmarks (+)")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "figures/template_scatter.png"), dpi=110)
        plt.close(fig)
        md += ["![template scatter](figures/template_scatter.png)", ""]
    else:
        md += ["Needs 3 or more subjects.", ""]
    if outl:
        md += ["Subject outliers: distance of a subject's consensus point to the mean of the other "
               "subjects' points, above both median + 3 MAD (of that distance over subjects) and 150 µm:", "",
               md_table(pd.DataFrame(outl)), ""]

    # 5. left/right symmetry
    md += ["## 7. Left/right symmetry", ""]
    pairs = {}
    for n in tpl:
        p = pair_of(n)
        if p:
            pairs.setdefault(p[0], {})[p[1]] = n
    pairs = {k: v for k, v in pairs.items() if len(v) == 2}
    if pairs:
        # template midsagittal plane: x = x0 (LPS x is left-right)
        x0 = np.median([(tpl[v["L"]][0] + tpl[v["R"]][0]) / 2 for v in pairs.values()])
        sym = []
        for base, v in sorted(pairs.items()):
            L, R = tpl[v["L"]], tpl[v["R"]]
            mirror = np.array([2 * x0 - R[0], R[1], R[2]])
            row = dict(pair=base.format("L/R"), template_asym_um=np.linalg.norm(L - mirror) * 1000)
            hL = cons[cons.landmark == v["L"]].set_index("scan")
            hR = cons[cons.landmark == v["R"]].set_index("scan")
            both = hL.index.intersection(hR.index)
            if len(both):
                cl = hL.loc[both, ["cx", "cy", "cz"]].to_numpy()
                cr = hR.loc[both, ["cx", "cy", "cz"]].to_numpy()
                cr_m = cr.copy()
                cr_m[:, 0] = 2 * x0 - cr[:, 0]
                row["subjects_asym_median_um"] = np.median(np.linalg.norm(cl - cr_m, axis=1)) * 1000
                row["consensus_err_L_um"] = hL.loc[both, "consensus_err_um"].median()
                row["consensus_err_R_um"] = hR.loc[both, "consensus_err_um"].median()
            sym.append(row)
        sym = pd.DataFrame(sym)
        sym.to_csv(os.path.join(out, "lr_symmetry.tsv"), sep="\t", index=False)
        sym_fig = None
        if "consensus_err_L_um" in sym:
            sp = sym.dropna(subset=["consensus_err_L_um"]).sort_values("consensus_err_L_um")
            fig, ax = plt.subplots(figsize=(7, 0.5 * len(sp) + 1.3))
            y = np.arange(len(sp))
            ax.hlines(y, sp["consensus_err_L_um"], sp["consensus_err_R_um"], color=GRID, lw=2)
            ax.plot(sp["consensus_err_L_um"], y, "o", color=C_L, ms=8, mec="white", mew=1.5, label="left")
            ax.plot(sp["consensus_err_R_um"], y, "o", color=C_R, ms=8, mec="white", mew=1.5, label="right")
            ax.set_yticks(y, sp["pair"])
            ax.set_xlabel("median consensus error over subjects (µm)")
            ax.set_title("Left vs right side of each landmark pair", loc="left")
            ax.grid(axis="y", visible=False)
            ax.legend(frameon=False, loc="lower right")
            fig.tight_layout()
            fig.savefig(os.path.join(out, "figures/lr_symmetry.png"), dpi=130)
            plt.close(fig)
            sym_fig = "![left/right](figures/lr_symmetry.png)"
        md += [
            f"Template midsagittal plane x = {x0:.3f} mm (LPS), from the paired template points. "
            "*Asym* = distance between the left point and the mirrored right point (0 for a perfectly "
            "symmetric pair). A large subject asymmetry, or one side much worse than the other, "
            "suggests an inconsistent placement on one side.",
            "",
            md_table(sym),
            "",
        ] + ([sym_fig, ""] if sym_fig else [])

    # 6. MRI vs lightsheet
    md += ["## 8. MRI vs lightsheet", ""]
    a = df[(df.path == "spim_to_template") & (df.inversion == "fixedpoint")]
    b = df[(df.path == "mri_to_template") & (df.inversion == "fixedpoint")]
    if not b.empty:
        k = ["scan", "method", "stain", "landmark"]
        m = a[k + ["tre_um"]].merge(b[k + ["tre_um"]], on=k, suffixes=("_ls", "_mri"))
        m["extra_from_mri_um"] = m["tre_um_mri"] - m["tre_um_ls"]
        t6 = m.groupby("landmark")[["tre_um_ls", "tre_um_mri", "extra_from_mri_um"]].median().sort_values(
            "extra_from_mri_um", ascending=False).reset_index()
        t6.to_csv(os.path.join(out, "mri_vs_ls.tsv"), sep="\t", index=False)
        ms = df[df.path == "mri_to_spim"].groupby("landmark")["tre_um"].median()
        t6["mri_to_spim_um"] = t6["landmark"].map(ms)
        fig, ax = plt.subplots(figsize=(6.2, 6))
        lim = [0, max(t6["tre_um_ls"].max(), t6["tre_um_mri"].max()) * 1.08]
        ax.plot(lim, lim, color=INK2, lw=1, ls="--")
        ax.text(lim[1] * 0.97, lim[1] * 0.93, "same error", ha="right", color=INK2, fontsize=8, rotation=45)
        ax.scatter(t6["tre_um_ls"], t6["tre_um_mri"], s=40, color=C_MRI, edgecolor="white", lw=1.2, zorder=3)
        placed = []  # label positions (data units), to step a label down past one it would hit
        for _, r in t6.nlargest(5, "extra_from_mri_um").iterrows():
            ly = r.tre_um_mri
            while any(abs(ly - y0) < lim[1] * 0.03 and abs(r.tre_um_ls - x0) < lim[1] * 0.3 for x0, y0 in placed):
                ly -= lim[1] * 0.03
            placed.append((r.tre_um_ls, ly))
            ax.annotate(r.landmark, (r.tre_um_ls, r.tre_um_mri), xytext=(r.tre_um_ls + lim[1] * 0.012, ly),
                        textcoords="data", fontsize=7, color=INK, va="center",
                        arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5) if ly != r.tre_um_mri else None)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_aspect("equal")
        ax.set_xlabel("TRE to template from the LS landmark (µm)")
        ax.set_ylabel("TRE to template from the MRI landmark (µm)")
        ax.set_title("Per landmark: MRI route vs LS route", loc="left")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "figures/mri_vs_ls.png"), dpi=130)
        plt.close(fig)
        md += [
            "Median TRE to the template from the lightsheet landmark vs from the MRI landmark (MRI→SPIM→"
            "template, same method and stain). A landmark whose MRI route is much worse than its "
            "lightsheet route, together with a high MRI→SPIM TRE, is worth checking on the MRI. "
            "In the figure each dot is a landmark; above the line the MRI route is worse "
            "(the largest gaps are named).",
            "",
            md_table(t6),
            "",
            "![MRI vs LS](figures/mri_vs_ls.png)",
            "",
        ]
    else:
        md += ["No MRI in this run.", ""]

    # 7. flags + screenshots
    flags = []
    for _, r in cons.iterrows():
        h = cons[cons.landmark == r.landmark]["consensus_err_um"]
        thr = max(np.median(h) + 3 * mad(h), 250) if len(h) >= 3 else 300
        if r.consensus_err_um > thr and r.ratio >= 2:
            flags.append(dict(scan=r.scan, subject=r.subject, landmark=r.landmark,
                              consensus_err_um=r.consensus_err_um, method_spread_um=r.method_spread_um,
                              reason="registrations agree with each other, not with the landmark"))
    for o in outl:
        flags.append(dict(scan="", subject=o["subject"], landmark=o["landmark"],
                          consensus_err_um=np.nan, method_spread_um=np.nan,
                          reason=f"outlier across subjects ({o['dist_to_others_um']:.0f} µm from the others)"))
    flags = pd.DataFrame(flags)
    md += ["## 9. Points worth a human look", ""]
    if flags.empty:
        md += ["None flagged.", ""]
    else:
        flags.to_csv(os.path.join(out, "flags.tsv"), sep="\t", index=False)
        md += [md_table(flags), ""]
        src = df[(df.path == PRIMARY[0]) & (df.inversion == PRIMARY[1])]
        for _, f in flags.drop_duplicates(["subject", "landmark"]).iterrows():
            s = src[(src.subject == f.subject) & (src.landmark == f.landmark)]
            c = cons[(cons.subject == f.subject) & (cons.landmark == f.landmark)]
            if s.empty or c.empty:
                continue
            p_src = xyz(s, "src")[0]
            imgs = sorted(glob.glob(os.path.join(run_root, f"sub-{f.subject}", "**", "micr",
                                                  "*level-5_desc-N4brain_SPIM.nii.gz"), recursive=True))
            fig, axs = plt.subplots(2, 3, figsize=(9, 6))
            if imgs:
                _slice_png(imgs[0], p_src, axs[0], f"SPIM ({os.path.basename(imgs[0]).split('stain-')[1].split('_')[0]})\n+ landmark")
            if tpl_image and os.path.exists(tpl_image):
                cpt = c[["cx", "cy", "cz"]].to_numpy()[0]
                _slice_png(tpl_image, tpl[f.landmark], axs[1], "template\n+ landmark, o consensus",
                           marks=[("cyan", cpt)])
            fig.suptitle(f"sub-{f.subject} {f.landmark}: {f.reason}", fontsize=9)
            fig.tight_layout()
            fn = f"figures/flag_sub-{f.subject}_{f.landmark}.png"
            fig.savefig(os.path.join(out, fn), dpi=100)
            plt.close(fig)
            md += [f"![{f.subject} {f.landmark}]({fn})", ""]

    # template-space markups: template points and every subject's consensus point
    mk = os.path.join(out, "markups")
    os.makedirs(mk, exist_ok=True)
    write_markups([(n, "template landmark", p) for n, p in sorted(tpl.items())],
                  os.path.join(mk, f"space-{template}_desc-template_landmarks.mrk.json"),
                  "template landmarks", [0.11, 0.69, 0.48])
    write_markups([(f"{r.subject}_{r.landmark}",
                    f"consensus of {r.n_regs} registrations, {r.consensus_err_um:.0f} um from the template point",
                    (r.cx, r.cy, r.cz)) for r in cons.itertuples()],
                  os.path.join(mk, f"space-{template}_desc-consensus_landmarks.mrk.json"),
                  "subject consensus", [0.92, 0.41, 0.20])
    md += ["## Files to view in 3D Slicer", "",
           f"`markups/space-{template}_desc-template_landmarks.mrk.json` (template points) and "
           f"`markups/space-{template}_desc-consensus_landmarks.mrk.json` (each subject's consensus point, "
           f"labelled `<subject>_<landmark>`), on the {template} template image. Every registration's warped "
           "points per scan are in `landmarks/sub-*/…/markups/` (listed in each `*_markups.tsv`).", ""]

    # summary figure, at the top of the report
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    ax = axs[0, 0]
    g = df[(df.path == PRIMARY[0]) & (df.inversion == PRIMARY[1])]
    ms = g.groupby("method")["tre_um"].agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)])
    ms.columns = ["med", "q25", "q75"]
    ms = ms.sort_values("med", ascending=False)
    y = np.arange(len(ms))
    ax.barh(y, ms["med"], height=0.55, color=C_LS)
    ax.errorbar(ms["med"], y, xerr=[ms["med"] - ms["q25"], ms["q75"] - ms["med"]], fmt="none", ecolor=INK2, lw=1)
    for yi, v in zip(y, ms["med"]):
        ax.text(6, yi, f"{v:.0f}", color="white", fontsize=8, va="center", fontweight="bold")
    ax.set_yticks(y, ms.index)
    ax.set_xlabel("median TRE (µm), all stains, bars = IQR")
    ax.set_title("A. LS → template by method", loc="left")
    ax.grid(axis="y", visible=False)
    if not failed.empty:
        fm = sorted(failed["method"].unique())
        some, never = [m for m in fm if m in ms.index], [m for m in fm if m not in ms.index]
        note = "; ".join(x for x in [
            ("did not converge on some scans (those left out): " + ", ".join(some)) if some else "",
            ("did not converge (left out): " + ", ".join(never)) if never else ""] if x)
        ax.text(0.99, 0.02, note, transform=ax.transAxes, ha="right", fontsize=7, color=INK2)
    ax = axs[0, 1]
    data = [best_rows_of(p, i)["tre_um"].to_numpy() for p, i, _ in paths]
    labs = [f"{lab}\n{best[(p, i)][0]}/{best[(p, i)][1]}" for p, i, lab in paths]
    bp = ax.boxplot(data, tick_labels=labs, widths=0.55, patch_artist=True, medianprops=dict(color=INK, lw=1.5),
                    flierprops=dict(ms=2, mec=INK2), whiskerprops=dict(color=INK2), capprops=dict(color=INK2))
    for b, (p, _, _) in zip(bp["boxes"], paths):
        col = C_MRI if "mri" in p else C_LS
        b.set(facecolor=col + "55", edgecolor=col)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_ylabel("TRE (µm)")
    ax.set_title("B. Each path, best method/stain", loc="left")
    ax.grid(axis="x", visible=False)
    ax = axs[1, 0]
    if subj_order:
        g = best_rows_of(*PRIMARY)
        sm_ = g.groupby("subject")["tre_um"].median().loc[subj_order]
        ax.plot(sm_.to_numpy(), np.arange(len(sm_)), "o", color=C_LS, ms=8, mec="white", mew=1.5, label="LS → template")
        if ("mri_to_template", "fixedpoint") in best:
            gm = best_rows_of("mri_to_template", "fixedpoint").groupby("subject")["tre_um"].median()
            gm = gm.reindex(subj_order)
            ax.plot(gm.to_numpy(), np.arange(len(gm)), "o", color=C_MRI, ms=8, mec="white", mew=1.5,
                    label="MRI → template")
            ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1, 1), ncol=2, fontsize=8)
        ax.set_yticks(np.arange(len(subj_order)), [f"sub-{s_}" for s_ in subj_order])
        ax.set_xlabel("median TRE (µm)")
        ax.set_xlim(left=0)
    ax.set_title("C. By subject", loc="left")
    ax.grid(axis="y", visible=False)
    ax = axs[1, 1]
    if region_order:
        g = best_rows_of(*PRIMARY)
        rm = g.groupby("region")["tre_um"].median().reindex(region_order[::-1])
        y = np.arange(len(rm))
        ax.barh(y, rm.to_numpy(), height=0.55, color=C_LS)
        for yi, v in zip(y, rm.to_numpy()):
            ax.text(v + 4, yi, f"{v:.0f}", va="center", color=INK, fontsize=8)
        ax.set_yticks(y, rm.index)
        ax.set_xlabel("median TRE (µm), LS → template")
    ax.set_title("D. By region", loc="left")
    ax.grid(axis="y", visible=False)
    nsubj = df["subject"].nunique()
    fig.suptitle(f"Landmark TRE summary: {nsubj} subjects, {df['landmark'].nunique()} landmarks",
                 x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "figures/summary.png"), dpi=130)
    plt.close(fig)
    md[header_end:header_end] = ["![summary](figures/summary.png)", ""]

    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(md))


if __name__ == "__main__":
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        ti = sm.input.template_image
        ti = ti if isinstance(ti, str) else (list(ti)[0] if len(ti) else None)
        main(sm.input.tsv, sm.params.template_landmarks, ti, sm.params.run_root,
             os.path.dirname(sm.output.report), sm.params.regions, sm.params.template)
    else:
        ap = argparse.ArgumentParser()
        ap.add_argument("--tsv", required=True)
        ap.add_argument("--template-landmarks", required=True)
        ap.add_argument("--template-image", default=None)
        ap.add_argument("--run-root", default=".")
        ap.add_argument("--out", required=True)
        ap.add_argument("--regions", default=None, help="TSV with columns region, pattern")
        ap.add_argument("--template", default="template")
        a = ap.parse_args()
        main(a.tsv, a.template_landmarks, a.template_image, a.run_root, a.out, a.regions, a.template)
