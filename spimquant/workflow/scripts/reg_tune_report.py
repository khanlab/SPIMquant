"""Side-by-side report comparing registration methods for one image pair."""

from regqc_lib import render_report

p = snakemake.params

candidates = [
    {
        "name": name,
        "warped": warped,
        "labels": labels,
        "warp": warp,
        "note": p.notes[name],
    }
    for name, warped, labels, warp in zip(
        p.methods,
        snakemake.input.warped,
        snakemake.input.labels,
        snakemake.input.warps,
    )
]

render_report(
    snakemake.output.report_html,
    title=f"{p.subject} to {p.template}",
    eyebrow="Registration method comparison",
    facts={"Moving": p.subject, "Fixed": p.template, "Stain": p.stain},
    fixed_path=snakemake.input.fixed,
    candidates=candidates,
    hint=f"Selected for this run: {p.selected}.",
)
