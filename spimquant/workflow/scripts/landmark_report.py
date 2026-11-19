"""Landmark report for --landmark_analysis: which tool/stain registers best, and
what the landmarks themselves look like. Flags points for a human look; never
changes them.

Run by the landmark_report rule, or:
  python landmark_report.py --tsv landmarks_tre_allscans.tsv --template-landmarks DIR \
      --template-image tpl-ABAv3_anat.nii.gz --run-root OUTPUT_DIR --out REPORT_DIR

Sections
 1. best method x stain per path (median TRE heatmaps)
 2. per-landmark TRE spread (boxplots)
 3. landmark vs registration: per subject x landmark, the consensus of all
    method x stain warps vs the landmark. If independent registrations agree with
    each other but not with the landmark, the landmark is the likely problem.
 4. consistency across subjects in template space (cluster spread, bias, outliers)
 5. left/right symmetry of paired landmarks
 6. MRI vs lightsheet agreement
 7. flagged points with screenshots
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


def main(tsv, tpl_dir, tpl_image, run_root, out):
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    df = pd.read_csv(tsv, sep="\t")
    df["scan"] = scan_id(df)
    df = df[df["inside_field"].astype(bool)].copy()
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

    # 2. per-landmark boxplots
    md += ["## 2. TRE per landmark", ""]
    for path, inv in [PRIMARY, ("mri_to_spim", "none"), ("mri_to_template", "fixedpoint")]:
        g = df[(df.path == path) & (df.inversion == inv)]
        if g.empty:
            continue
        order = g.groupby("landmark")["tre_um"].median().sort_values().index
        fig, ax = plt.subplots(figsize=(10, 0.28 * len(order) + 1.2))
        ax.boxplot([g[g.landmark == n]["tre_um"] for n in order], vert=False, labels=order,
                   flierprops=dict(ms=2))
        ax.set_xlabel("TRE (µm), all subjects x methods x stains")
        ax.set_title(f"{path} [{inv}]", fontsize=9)
        ax.tick_params(axis="y", labelsize=7)
        fig.tight_layout()
        fn = f"figures/box_{path}_{inv}.png"
        fig.savefig(os.path.join(out, fn), dpi=110)
        plt.close(fig)
        md += [f"![{path}]({fn})", ""]

    # 3. landmark vs registration (consensus)
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
        "## 3. Landmark or registration?",
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
    md += ["## 4. Consistency across subjects (template space)", ""]
    rows4, outl = [], []
    if cons["subject"].nunique() >= 3:
        for name, h in cons.groupby("landmark"):
            C = h[["cx", "cy", "cz"]].to_numpy()
            t = h[["tx", "ty", "tz"]].to_numpy()[0]
            m = C.mean(0)
            d = np.linalg.norm(C - m, axis=1) * 1000
            rows4.append(dict(landmark=name, n=len(h), spread_rms_um=np.sqrt((d ** 2).mean()),
                              bias_of_mean_um=np.linalg.norm(m - t) * 1000))
            thr = max(np.median(d) + 3 * mad(d), 150)
            for (_, r), di in zip(h.iterrows(), d):
                if di > thr:
                    outl.append(dict(subject=r.subject, landmark=name, dist_to_cluster_um=di, threshold_um=thr))
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
        md += ["Subject outliers (distance to the cluster of the other subjects > median + 3 MAD and > 150 µm):", "",
               md_table(pd.DataFrame(outl)), ""]

    # 5. left/right symmetry
    md += ["## 5. Left/right symmetry", ""]
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
        md += [
            f"Template midsagittal plane x = {x0:.3f} mm (LPS), from the paired template points. "
            "*Asym* = distance between the left point and the mirrored right point (0 for a perfectly "
            "symmetric pair). A large subject asymmetry, or one side much worse than the other, "
            "suggests an inconsistent placement on one side.",
            "",
            md_table(sym),
            "",
        ]

    # 6. MRI vs lightsheet
    md += ["## 6. MRI vs lightsheet", ""]
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
        md += [
            "Median TRE to the template from the lightsheet landmark vs from the MRI landmark (MRI→SPIM→"
            "template, same method and stain). A landmark whose MRI route is much worse than its "
            "lightsheet route, together with a high MRI→SPIM TRE, is worth checking on the MRI.",
            "",
            md_table(t6),
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
                          reason=f"outlier across subjects ({o['dist_to_cluster_um']:.0f} µm from cluster)"))
    flags = pd.DataFrame(flags)
    md += ["## 7. Points worth a human look", ""]
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

    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(md))


if __name__ == "__main__":
    if "snakemake" in globals():
        sm = globals()["snakemake"]
        ti = sm.input.template_image
        ti = ti if isinstance(ti, str) else (list(ti)[0] if len(ti) else None)
        main(sm.input.tsv, sm.params.template_landmarks, ti, sm.params.run_root,
             os.path.dirname(sm.output.report))
    else:
        ap = argparse.ArgumentParser()
        ap.add_argument("--tsv", required=True)
        ap.add_argument("--template-landmarks", required=True)
        ap.add_argument("--template-image", default=None)
        ap.add_argument("--run-root", default=".")
        ap.add_argument("--out", required=True)
        a = ap.parse_args()
        main(a.tsv, a.template_landmarks, a.template_image, a.run_root, a.out)
