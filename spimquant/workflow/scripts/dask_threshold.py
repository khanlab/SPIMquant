import numpy as np
import json
from zarrnii import ZarrNii
from lib.utils import get_zarr_store

if snakemake.config['use_coiled']:
    from coiled import Cluster
    cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=30,idle_timeout='1 hour')
    client = cluster.get_client()


store = get_zarr_store(snakemake.params.spim_n4_uri)

in_orient = snakemake.config['in_orientation']
orient_opt = {} if in_orient == None else {'orientation': in_orient}


#we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
# but we need to adjust the z_level_offset if the level is not zero -- TODO: deprecate z_level_offset in zarrnii by automatically calculating it based on the xy to z ratio.. 
z_level_offset=min(0,-2+int(snakemake.wildcards.level))

znimg_hires = ZarrNii.from_ome_zarr(store,z_level_offset=z_level_offset,**orient_opt)

def threshold_block(x):
    return np.where(x > snakemake.params.threshold,100,0) 


#now, we perform thresholding on hires, and save the result in a new ome-zarr 
znimg_hires.darr = znimg_hires.darr.map_blocks(threshold_block,dtype=np.uint8,meta=np.array((), dtype=np.uint8))
znimg_hires.to_ome_zarr(snakemake.params.mask_uri,max_layer=5)

