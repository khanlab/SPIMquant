from coiled import Cluster
import zarr
import dask.array as da
import numpy as np
from skimage.filters import threshold_multiotsu
from zarrnii import ZarrNii

cluster = Cluster(n_workers=[4,20],spot_policy='spot_with_fallback')
client = cluster.get_client()


hires_level=int(snakemake.wildcards.level)
ds_level=int(snakemake.wildcards.dslevel)
ds_level_z=ds_level-1 #z downsampling one less (since already lower-res


# use downsampled level to get globally optimum threshold
znimg_ds = ZarrNii.from_path(snakemake.params.spim_n4_uri,level=ds_level).downsample(along_z=2**ds_level_z)
znimg_n4 = ZarrNii.from_path(snakemake.params.spim_n4_uri,level=hires_level)

print('n4 corr')
print(znimg_n4.darr.shape)
print(znimg_n4.darr.chunks)
print(znimg_n4.darr.blocks.shape)



multi_thresholds = threshold_multiotsu(znimg_ds.darr.compute(),classes=snakemake.params.otsu_n_classes)

def threshold_block(x):
    """ our thresholding function, returns 100 if above, 0 if not"""
    return np.where(x > multi_thresholds[snakemake.params.otsu_threshold_index], 100, 0)



#now, we perform thresholding on hires, and save the result in a new ome-zarr 
znimg_n4.darr = znimg_n4.darr.map_blocks(threshold_block,dtype=np.uint8,meta=np.array((), dtype=np.uint8))

znimg_n4.to_ome_zarr(snakemake.params.mask_uri,max_layer=5)

