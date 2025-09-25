from dask.diagnostics import ProgressBar
from lib.utils import get_zarr_store

from zarrnii import ZarrNii

store = get_zarr_store(snakemake.params.uri)

znimg = ZarrNii.from_ome_zarr(
    store,
    level=int(snakemake.wildcards.level),
    channel_labels=[snakemake.wildcards.stain],
    downsample_near_isotropic=True,
)

with ProgressBar():
    znimg.to_nifti(snakemake.output.nii)
