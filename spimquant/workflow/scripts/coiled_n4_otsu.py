from coiled import Cluster
import zarr
import dask.array as da
import dask
from dask.array.overlap import overlap, trim_overlap
import numpy as np
from math import sqrt
from skimage.feature import blob_dog
from skimage.filters import threshold_multiotsu
import sparse
import coiled
from zarrnii import ZarrNii

cluster = Cluster(name=snakemake.params.cluster_name,n_workers=[4,20],spot_policy='spot_with_fallback')
client = cluster.get_client()


hires_level=snakemake.wildcards.level
ds_level=snakemake.wildcards.dslevel
ds_level_z=ds_level-1 #z downsampling one less (since already lower-res
stain=snakemake.wildcards.stain



zi = zarr.open(snakemake.params.spim_uri)
attrs=zi['/'].attrs.asdict()

#get channel index from omero metadata
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(stain)


# # apply N4 bias field to full-res (level 0) image


# ### TODO: figure out optimal chunk size for this ahead of time:
znimg_hires = ZarrNii.from_path(snakemake.params.spim_uri,channels=[channel_index],level=hires_level,chunks=(1,160,320,320),rechunk=True) #chunk size comes from upsampled array




#ok, now we have the bias field.. let's resample it to the level 
# where we perform thresholding

#first, we write the nifti n4 bias field to ome zarr
ZarrNii.from_path(snakemake.input.n4_bf_ds,chunks=(1,20,20,20),rechunk=True).to_ome_zarr(snakemake.params.bf_ds_uri)


#then we upsample the bias field (TODO: could apply downsampled bias field in map_blocks instead (using zoom in there) would be one less step)
znimg_biasfield_upsampled = ZarrNii.from_path(snakemake.params.bf_ds_uri).upsample(to_shape=znimg_hires.darr.shape)
znimg_biasfield_upsampled.to_ome_zarr(snakemake.params.bf_us_uri,max_layers=0)


#now multiply biasfield and image together, 
znimg_hires.darr = znimg_hires.darr / ZarrNii.from_path(snakemake.params.bf_us_uri).darr


#write the full-res N4-corrected OME-zarr dataset
znimg_hires.to_ome_zarr(snakemake.params.spim_n4_uri,max_layers=5)

    
#TODO -- consider breaking up script here, separateing N4 correction from thresholding.. maybe later when we are optimizing N4 and seg params..

# use downsampled level to get globally optimum threshold
znimg_ds = ZarrNii.from_path(snakemake.params.spim_n4_uri,level=ds_level).downsample(along_z=2**ds_level_z)
znimg_n4 = ZarrNii.from_path(snakemake.params.spim_n4_uri,level=hires_level)


multi_thresholds = threshold_multiotsu(znimg_ds.darr.compute(),classes=snakemake.params.otsu_n_classes)

def threshold_block(x):
    """ our thresholding function, returns 100 if above, 0 if not"""
    return np.where(x > multi_thresholds[snakemake.params.otsu_threshold_index], 100, 0)



#now, we perform thresholding on hires, and save the result in a new ome-zarr 
znimg_n4.darr = znimg_n4.darr.map_blocks(threshold_block,dtype=np.uint8,meta=np.array((), dtype=np.uint8))
znimg_n4.to_ome_zarr(snakemake.params.mask_uri,max_layer=5)


