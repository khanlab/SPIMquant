"""
Generic rule to create BIDS JSON sidecar files for tabular outputs.

One rule covers both TSV and Parquet files. The JSON sidecar is written
alongside the tabular file with the same stem and a ``.json`` extension.
"""


rule create_tabular_json_sidecar:
    """Create a BIDS JSON sidecar describing the columns of a tabular output."""
    wildcard_constraints:
        tabular_ext="tsv|parquet",
    input:
        "{prefix}.{tabular_ext}",
    output:
        "{prefix}.json",
    params:
        column_descriptions=config["column_descriptions"],
        stats_maps=config["stats_maps"],
    threads: 1
    resources:
        mem_mb=500,
        runtime=5,
    script:
        "../scripts/create_json_sidecar.py"
