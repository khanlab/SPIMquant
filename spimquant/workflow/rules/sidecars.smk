"""
Generic rule to create BIDS JSON sidecar files for tabular outputs.

One rule covers both TSV and Parquet files. The JSON sidecar is written
alongside the tabular file with the same stem and a ``.json`` extension.
"""


rule create_tabular_json_sidecar_from_tsv:
    """Create a BIDS JSON sidecar describing the columns of a tabular output."""
    input:
        "{prefix}.tsv",
    output:
        "{prefix}.json",
    threads: 1
    resources:
        mem_mb=500,
        runtime=5,
    params:
        column_descriptions=config["column_descriptions"],
        stats_maps=config["stats_maps"],
    script:
        "../scripts/create_json_sidecar.py"


rule create_tabular_json_sidecar_from_parquet:
    """Create a BIDS JSON sidecar describing the columns of a tabular output."""
    input:
        "{prefix}.parquet",
    output:
        "{prefix}.json",
    threads: 1
    resources:
        mem_mb=500,
        runtime=5,
    params:
        column_descriptions=config["column_descriptions"],
        stats_maps=config["stats_maps"],
    script:
        "../scripts/create_json_sidecar.py"
