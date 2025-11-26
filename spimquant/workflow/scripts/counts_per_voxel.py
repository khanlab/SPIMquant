import numpy as np
from zarrnii import ZarrNii, density_from_points
from dask.diagnostics import ProgressBar
import pandas as pd

level = int(snakemake.wildcards.level)
stain = snakemake.wildcards.stain

img = ZarrNii.from_ome_zarr(
    snakemake.input.ref_spim,
    level=level,
    channel_labels=[stain],
    downsample_near_isotropic=True,
)


df = pd.read_parquet(snakemake.input.regionprops_parquet)

points = df[snakemake.params.coord_column_names].values

# Create counts map (zarrnii is calling this density right now)..
counts = density_from_points(points, img, in_physical_space=True)
with ProgressBar():
    counts.to_nifti(snakemake.output.counts_nii)
