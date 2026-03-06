"""Map mask colocalization results to atlas regions.

This script processes the mask colocalization parquet (instance segmentation objects
annotated with in_mask=0/1) and maps the objects that are inside the mask (in_mask=1)
to anatomical atlas regions using template-space coordinates.

The script generates two output TSV files:
- maskcoloc_tsv: Contains per-region statistics of objects inside the mask, including
  region labels and colocalization metrics (counts of in-mask objects per region).
- counts_tsv: Contains total count statistics for each atlas region (all instance
  objects, not filtered by in_mask).

This script is intended to be run as part of the SPIMquant Snakemake workflow, and
expects the `snakemake` object to be available for accessing input/output file paths
and parameters.

Inputs (via snakemake.input):
    - maskcoloc_parquet: Parquet file with annotated instance regionprops and in_mask column.
    - dseg: Atlas segmentation file (in template space).
    - label_tsv: TSV file with atlas region labels.

Outputs (via snakemake.output):
    - maskcoloc_tsv: TSV file with per-region mask colocalization statistics
      (only objects inside mask, in_mask=1).
    - counts_tsv: TSV file with per-region counts of all instance objects.

Parameters (via snakemake.params):
    - coord_column_names: List of column names for template-space coordinates
      (e.g. ['template_x', 'template_y', 'template_z']).
"""

import pandas as pd
from zarrnii import ZarrNiiAtlas

# Load the mask colocalization results
df_full = pd.read_parquet(snakemake.input.maskcoloc_parquet)

# Load the atlas
atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

mask_stain = df_full["mask_stain"].iloc[0] if len(df_full) > 0 else "unknown"
n_in_mask = int(df_full["in_mask"].sum()) if len(df_full) > 0 else 0

if len(df_full) > 0:
    # Map ALL instance objects to atlas regions (for the counts output)
    all_objects = df_full.to_dict(orient="list")
    _, df_counts = atlas.label_region_properties(
        all_objects,
        coord_column_names=snakemake.params.coord_column_names,
        include_names=True,
    )
else:
    # No instance objects: create an empty counts DataFrame
    df_counts = pd.DataFrame()

# Filter to objects inside the mask (in_mask=1) for the colocalization stats
df_in_mask = df_full[df_full["in_mask"] == 1] if len(df_full) > 0 else df_full

if len(df_in_mask) > 0:
    in_mask_objects = df_in_mask.to_dict(orient="list")
    df_maskcoloc, _ = atlas.label_region_properties(
        in_mask_objects,
        coord_column_names=snakemake.params.coord_column_names,
        include_names=True,
    )
else:
    # No objects inside the mask: output an empty DataFrame.
    # Use df_counts as the basis if it is non-empty, otherwise create an empty one.
    df_maskcoloc = df_counts.iloc[0:0].copy() if len(df_counts) > 0 else pd.DataFrame()

df_maskcoloc.to_csv(snakemake.output.maskcoloc_tsv, sep="\t", index=False)
df_counts.to_csv(snakemake.output.counts_tsv, sep="\t", index=False)

print("Mask colocalization atlas mapping complete:")
print(f"  Mask stain: {mask_stain}")
print(f"  Total instance objects mapped: {len(df_full)}")
print(f"  Objects inside mask mapped: {n_in_mask}")
