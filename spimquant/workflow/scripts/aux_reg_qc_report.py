"""QC report for the registration of an additional scan to its reference scan."""

import json

from regqc_lib import render_report

p = snakemake.params
inp = snakemake.input

with open(inp.init_info) as f:
    init = json.load(f)

final_name = {"rigid": "Rigid", "affine": "Affine", "deform": "Deformable"}[p.method]
init_note = (
    f"match, NCC {init['ncc']:.3f}" if init.get("ncc") is not None else init["method"]
)

render_report(
    snakemake.output.report_html,
    title=f"{p.ref} to {p.scan}",
    eyebrow="Additional scan registration QC",
    facts={
        "Moving (reference scan)": p.ref,
        "Fixed (additional scan)": p.scan,
        "Stains": ", ".join(p.stains),
        "Initialization": init_note,
        "Initial translation (mm)": ", ".join(
            f"{v:+.2f}" for v in init["translation_mm"]
        ),
        "Method": p.method,
        "Parameters": p.parameters,
    },
    fixed_path=inp.fixed[0],
    candidates=[
        {"name": "Initial", "warped": inp.warped_init[0], "note": init_note},
        {
            "name": final_name,
            "warped": inp.warped_final[0],
            "labels": inp.dseg,
            "warp": inp.get("warp"),
        },
    ],
    hint=(
        "Initial is the reference scan placed by the initial translation alone "
        "(--aux_reg_init): a poor match there (low NCC) can leave the registration "
        "in the wrong place. Atlas outlines are the labels this scan gets through "
        "the reference scan."
    ),
    panel_width=620,
    stains=[
        {"name": name, "fixed": fixed, "warped": [w_init, w_final]}
        for name, fixed, w_init, w_final in zip(
            p.stains, inp.fixed, inp.warped_init, inp.warped_final
        )
    ],
)
