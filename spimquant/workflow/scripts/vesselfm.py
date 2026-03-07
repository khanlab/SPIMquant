from zarrnii import ZarrNii
from vesselfm.zarrnii_plugin import VesselFMPlugin

import dask
dask.config.set(scheduler="threads", num_workers=snakemake.threads)

from dask.diagnostics import ProgressBar

znimg = ZarrNii.from_ome_zarr(
    snakemake.input.spim,
    level=int(snakemake.wildcards.level),
    channel_labels=[snakemake.wildcards.stain],
)
znimg_mask = znimg.segment(VesselFMPlugin, chunk_size=(1, 128, 128, 128))

znimg_mask = znimg_mask * 100

with ProgressBar():
    znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5, zarr_format=2)
