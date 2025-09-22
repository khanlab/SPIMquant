import numpy as np
import json
from zarrnii import ZarrNii
from lib.utils import get_zarr_store

if snakemake.config["use_coiled"]:
    from coiled import Cluster

    cluster = Cluster(
        name="coiled-snakemake",
        package_sync_ignore=["spimquant"],
        n_workers=30,
        idle_timeout="1 hour",
    )
    client = cluster.get_client()


store = get_zarr_store(snakemake.params.spim_n4_uri)

znimg_hires = ZarrNii.from_ome_zarr(store)


def threshold_block(x):
    return np.where(x > snakemake.params.threshold, 100, 0)


# now, we perform thresholding on hires, and save the result in a new ome-zarr
znimg_hires.darr = znimg_hires.darr.map_blocks(
    threshold_block, dtype=np.uint8, meta=np.array((), dtype=np.uint8)
)
znimg_hires.to_ome_zarr(snakemake.params.mask_uri, max_layer=5)
