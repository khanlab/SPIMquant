import numpy as np
from zarrnii import ZarrNii, density_from_points
from dask.diagnostics import ProgressBar

level = int(snakemake.wildcards.dslevel)
stain = snakemake.wildcards.stain

img = ZarrNii.from_ome_zarr(
    snakemake.input.ref_spim,
    level=level,
    channel_labels=[stain],
    downsample_near_isotropic=True,
)

# Create counts map (zarrnii is calling this density right now)..
counts = density_from_points(
    snakemake.input.centroids_parquet, img, in_physical_space=True
)
with ProgressBar():
    counts.to_nifti(snakemake.output.counts_nii)
