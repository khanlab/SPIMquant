from coiled import Cluster
import zarr
import dask.array as da
import numpy as np
from zarrnii import ZarrNii

cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=30,idle_timeout='1 hour')
client = cluster.get_client()


hires_level=int(snakemake.wildcards.level)


# use downsampled level to get globally optimum threshold
znimg_hires = ZarrNii.from_ome_zarr(snakemake.params.spim_n4_uri,level=hires_level)

print(znimg_hires.darr.shape)
print(znimg_hires.darr.chunks)
print(znimg_hires.darr.blocks.shape)


multi_thresholds = np.load(snakemake.input.otsu_thresholds)


def threshold_block(x):
    """ our thresholding function, returns 100 if above, 0 if not"""
    return np.where(x > multi_thresholds[snakemake.params.otsu_threshold_index], 100, 0)



#now, we perform thresholding on hires, and save the result in a new ome-zarr 
znimg_hires.darr = znimg_hires.darr.map_blocks(threshold_block,dtype=np.uint8,meta=np.array((), dtype=np.uint8))

znimg_hires.to_ome_zarr(snakemake.params.mask_uri,max_layer=5)

