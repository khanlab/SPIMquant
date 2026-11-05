"""Side-by-side report comparing registration methods for one image pair."""

from regqc_lib import render_report

p = snakemake.params
inp = snakemake.input
n_methods = len(p.methods)
warped = list(inp.warped)
fixed = [inp.fixed, *inp.extra_fixed]

# inp.warped has an image per method for each stain in turn, or, if the moving image
# has one contrast, an image per method that is compared to every stain
per_stain = len(warped) == n_methods * len(fixed)
stains = [
    {
        "name": name,
        "fixed": fixed_path,
        "warped": warped[s * n_methods : (s + 1) * n_methods]
        if per_stain
        else warped[:n_methods],
    }
    for s, (name, fixed_path) in enumerate(zip(p.stains, fixed))
]
candidates = [
    {
        "name": name,
        "warped": warped[c],
        "labels": labels,
        "warp": warp,
        "note": p.notes[name],
    }
    for c, (name, labels, warp) in enumerate(zip(p.methods, inp.labels, inp.warps))
]

render_report(
    snakemake.output.report_html,
    title=f"{p.subject} to {p.template}",
    eyebrow="Registration method comparison",
    facts={"Moving": p.subject, "Fixed": p.template, "Stains": ", ".join(p.stains)},
    fixed_path=inp.fixed,
    candidates=candidates,
    hint=f"Selected for this run: {p.selected}.",
    default_layer=getattr(p, "default_layer", "warped"),
    stains=stains,
)
