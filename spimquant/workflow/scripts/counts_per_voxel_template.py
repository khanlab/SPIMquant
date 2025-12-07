import numpy as np
from zarrnii import ZarrNii, density_from_points
from dask.diagnostics import ProgressBar
import pandas as pd

stain = snakemake.wildcards.stain

img = ZarrNii.from_nifti(
    snakemake.input.template,
)


df = pd.read_parquet(snakemake.input.regionprops_parquet)

df = df[df['stain'] == snakemake.wildcards.stain]
points = df[snakemake.params.coord_column_names].values

# Create counts map (zarrnii is calling this density right now)..
counts = density_from_points(points, img, in_physical_space=True)
with ProgressBar():
    counts.to_nifti(snakemake.output.counts_nii)
