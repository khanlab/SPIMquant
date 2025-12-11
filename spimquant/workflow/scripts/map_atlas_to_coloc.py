"""Map computed colocalized objects to atlas regions and generate labeled statistics.

This script processes colocalization results (a parquet file containing colocalized object pairs
with spatial coordinates) and maps these objects to anatomical atlas regions using template-space
coordinates. The mapping is performed using a ZarrNiiAtlas object, which reads the atlas segmentation
and label information.

The script generates two output TSV files:
- coloc_tsv: Contains per-region statistics of colocalized objects, including region labels and counts.
- counts_tsv: Contains count statistics for each atlas region, regardless of colocalization.

This script is intended to be run as part of the SPIMquant Snakemake workflow, and expects the
`snakemake` object to be available for accessing input and output file paths, as well as parameters.

Inputs (via snakemake.input):
- coloc_parquet: Parquet file with colocalized object pairs and their coordinates.
- dseg: Atlas segmentation file (in template space).
- label_tsv: TSV file with atlas region labels.

Outputs (via snakemake.output):
- coloc_tsv: TSV file with per-region colocalization statistics.
- counts_tsv: TSV file with per-region object counts.

Parameters (via snakemake.params):
- coord_column_names: List of column names specifying the coordinate columns in the colocalization data.
"""

from zarrnii import ZarrNiiAtlas
import pandas as pd

coloc = pd.read_parquet(snakemake.input.coloc_parquet).to_dict(orient="list")
atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

df_coloc, df_counts = atlas.label_region_properties(
    coloc,
    coord_column_names=snakemake.params.coord_column_names,
    include_names=True,
)

df_coloc.to_csv(snakemake.output.coloc_tsv, sep="\t", index=False)
df_counts.to_csv(snakemake.output.counts_tsv, sep="\t", index=False)
