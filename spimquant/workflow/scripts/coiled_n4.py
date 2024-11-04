from coiled import Cluster
import zarr
import dask.array as da
from zarrnii import ZarrNii

cluster = Cluster(package_sync_ignore=['spimquant'],n_workers=[4,20],spot_policy='spot_with_fallback')
client = cluster.get_client()


hires_level=int(snakemake.wildcards.level)
ds_level=int(snakemake.wildcards.dslevel)

stain=snakemake.wildcards.stain

zi = zarr.open(snakemake.params.spim_uri)
attrs=zi['/'].attrs.asdict()

#get channel index from omero metadata
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(stain)

# apply N4 bias field to full-res (level 0) image

znimg_hires = ZarrNii.from_path(snakemake.params.spim_uri,channels=[channel_index],level=hires_level)

print(znimg_hires.darr)
#ok, now we have the bias field.. let's resample it to the level 
# where we perform thresholding

#first, we write the nifti n4 bias field to ome zarr
ZarrNii.from_path(snakemake.input.n4_bf_ds,chunks=(1,8,8,4),rechunk=True).to_ome_zarr(snakemake.params.bf_ds_uri) #chunks here are nifti, so c,x,y,z
znimg_biasfield_downsampled =  ZarrNii.from_path(snakemake.params.bf_ds_uri)
print(znimg_biasfield_downsampled.darr)

znimg_n4 = znimg_hires.divide_by_downsampled(znimg_biasfield_downsampled)


print(znimg_n4.darr)

#write the full-res N4-corrected OME-zarr dataset
znimg_n4.to_ome_zarr(snakemake.params.spim_n4_uri,max_layers=6)
    

