"""QC report for one registration: its affine and deformable stages."""

from regqc_lib import render_report

p = snakemake.params
# labels warped with the deformable registration do not apply to the affine stage
labels_warped = getattr(p, "labels_warped", False)

render_report(
    snakemake.output.report_html,
    title=f"{p.subject} to {p.template}",
    eyebrow="Registration QC",
    facts={
        "Moving": p.subject,
        "Fixed": p.template,
        "Stain": p.stain,
        "Method": p.method,
        "Parameters": p.parameters,
    },
    fixed_path=snakemake.input.template,
    candidates=[
        {
            "name": "Affine",
            "warped": snakemake.input.warped_affine,
            "labels": None if labels_warped else snakemake.input.dseg,
        },
        {
            "name": "Deformable",
            "warped": snakemake.input.warped_deform,
            "labels": snakemake.input.dseg,
            "warp": snakemake.input.warp,
        },
    ],
    panel_width=620,
    default_layer="fixed" if labels_warped else "warped",
)
