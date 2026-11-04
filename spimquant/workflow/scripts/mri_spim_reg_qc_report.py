"""QC report for the registration of a subject's MRI to its SPIM image."""

from regqc_lib import render_report

p = snakemake.params

render_report(
    snakemake.output.report_html,
    title=f"MRI to SPIM, sub-{snakemake.wildcards.subject}",
    eyebrow="Registration QC",
    facts={"Moving": "MRI", "Fixed": "SPIM", "Stain": p.stain_for_reg},
    fixed_path=snakemake.input.spim,
    candidates=[
        {"name": "Linear", "warped": snakemake.input.warped_affine},
        {
            "name": "Deformable",
            "warped": snakemake.input.warped_deform,
            "warp": snakemake.input.warp,
        },
    ],
    panel_width=620,
)
