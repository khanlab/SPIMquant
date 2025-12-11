import numpy as np
from zarrnii import ZarrNii, density_from_points
from dask.diagnostics import ProgressBar
import dask
import pandas as pd

img = ZarrNii.from_nifti(
    snakemake.input.template,
)

if hasattr(snakemake.wildcards, "level"):
    img = img.downsample(level=int(snakemake.wildcards.level))

dask.config.set(scheduler="threads", num_workers=snakemake.threads)

df = pd.read_parquet(snakemake.input.coloc_parquet)

points = df[snakemake.params.coord_column_names].values

# Create counts map (zarrnii is calling this density right now)..
counts = density_from_points(points, img, in_physical_space=True)
with ProgressBar():
    counts.to_nifti(snakemake.output.counts_nii)
