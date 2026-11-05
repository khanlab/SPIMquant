"""QC report for one registration: its affine and deformable stages."""

from regqc_lib import render_report

p = snakemake.params
inp = snakemake.input
# labels warped with the deformable registration do not apply to the affine stage
labels_warped = getattr(p, "labels_warped", False)

# one row per stain. Without extra_warped, the moving image has a single contrast,
# which is compared to each stain; with it, subject stains are warped separately,
# and only the stain the registration was checked on has an affine stage to show
extra_warped = list(inp.get("extra_warped", []))
stains = [
    {
        "name": p.stains[0],
        "fixed": inp.template,
        "warped": [inp.warped_affine, inp.warped_deform],
    }
]
for i, (name, fixed) in enumerate(zip(p.stains[1:], inp.extra_fixed)):
    stains.append(
        {
            "name": name,
            "fixed": fixed,
            "warped": [None, extra_warped[i]]
            if extra_warped
            else [inp.warped_affine, inp.warped_deform],
        }
    )

render_report(
    snakemake.output.report_html,
    title=f"{p.subject} to {p.template}",
    eyebrow="Registration QC",
    facts={
        "Moving": p.subject,
        "Fixed": p.template,
        "Stains": ", ".join(p.stains),
        "Method": p.method,
        "Parameters": p.parameters,
    },
    fixed_path=inp.template,
    candidates=[
        {
            "name": "Affine",
            "warped": inp.warped_affine,
            "labels": None if labels_warped else inp.dseg,
        },
        {
            "name": "Deformable",
            "warped": inp.warped_deform,
            "labels": inp.dseg,
            "warp": inp.warp,
        },
    ],
    hint=(
        f"The affine stage is only shown for {p.stains[0]}." if extra_warped else ""
    ),
    panel_width=620,
    default_layer="fixed" if labels_warped else "warped",
    stains=stains,
)
