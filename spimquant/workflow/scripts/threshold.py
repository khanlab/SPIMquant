from dask_setup import get_dask_client
from zarrnii import ZarrNii

with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

    znimg_hires = ZarrNii.from_ome_zarr(
        snakemake.input.corrected, **snakemake.params.zarrnii_kwargs
    )

    print("thresholding image, saving as ome zarr")
    znimg_mask = znimg_hires.segment_threshold(snakemake.params.threshold)

    # multiplying binary mask by 100 (so values are 0  and 100) to enable
    # field fraction calculation by subsequent local-mean downsampling
    znimg_mask = znimg_mask * 100

    # write to ome_zarr
    znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5)
