"""One HTML report on all scans of a dataset, from their qc_scan.py outputs."""

import json

import pandas as pd
from qcreport_lib import build_report_data, metrics_table, render_report

scans = []
for path in snakemake.input.scans:
    with open(path) as f:
        scans.append(json.load(f))

participants = None
if snakemake.input.participants:
    participants = pd.read_csv(
        str(snakemake.input.participants), sep="\t", dtype=str
    ).fillna("")

data = build_report_data(
    scans,
    title=snakemake.params.title,
    stain_aliases=snakemake.params.stain_aliases,
    participants=participants,
    group_by=snakemake.params.group_by,
)
metrics_table(data).to_csv(snakemake.output.tsv, sep="\t", index=False)
render_report(snakemake.output.html, data)
