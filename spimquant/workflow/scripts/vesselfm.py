from zarrnii import ZarrNii
from vesselfm.zarrnii_plugin import VesselFMPlugin

from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client

with get_dask_client("threads", snakemake.threads):
    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.spim,
        level=int(snakemake.wildcards.level),
        channel_labels=[snakemake.wildcards.stain],
        downsample_near_isotropic=True,
        **snakemake.params.zarrnii_kwargs,
    )
    znimg_mask = znimg.segment(VesselFMPlugin, **snakemake.params.vesselfm_kwargs)

    znimg_mask = znimg_mask * 100

    with ProgressBar():
        znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5, zarr_format=2)
