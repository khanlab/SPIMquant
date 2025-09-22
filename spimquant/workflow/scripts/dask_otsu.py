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

# we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
znimg_hires = ZarrNii.from_ome_zarr(store)

with open(snakemake.input.otsu_thresholds, "r") as f:
    multi_thresholds = json.load(f)


def threshold_block(x):
    """use the thresholds to discretize, which np.digitize does nicely, then we set the chosen
    label index to 100 and others to 0"""
    thresholds = multi_thresholds[str(snakemake.params.otsu_k)]
    return np.where(
        np.digitize(x, thresholds) == snakemake.params.otsu_threshold_index, 100, 0
    )


# now, we perform thresholding on hires, and save the result in a new ome-zarr
znimg_hires.darr = znimg_hires.darr.map_blocks(
    threshold_block, dtype=np.uint8, meta=np.array((), dtype=np.uint8)
)
znimg_hires.to_ome_zarr(snakemake.params.mask_uri, max_layer=5)
