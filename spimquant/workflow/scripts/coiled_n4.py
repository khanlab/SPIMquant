from coiled import Cluster
import zarr
import dask.array as da
from zarrnii import ZarrNii


hires_level=int(snakemake.wildcards.level)
ds_level=int(snakemake.wildcards.dslevel)

stain=snakemake.wildcards.stain

zi = zarr.open(snakemake.params.spim_uri)
attrs=zi['/'].attrs.asdict()

#get channel index from omero metadata
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(stain)



# # apply N4 bias field to full-res (level 0) image

hires_shape = ZarrNii.from_path(snakemake.params.spim_uri,channels=[channel_index],level=hires_level).darr.shape

#ok, now we have the bias field.. let's resample it to the level
# where we perform thresholding

#first, we write the nifti n4 bias field to ome zarr
ZarrNii.from_path(snakemake.input.n4_bf_ds,chunks=(1,25,25,25),rechunk=True).to_ome_zarr(snakemake.params.bf_ds_uri)


cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=[4,20])
client = cluster.get_client()



#then we upsample the bias field
znimg_biasfield_upsampled = ZarrNii.from_path(snakemake.params.bf_ds_uri).upsample(to_shape=hires_shape)

print('biasfield_upsampled')
print(znimg_biasfield_upsampled.darr.shape)
print(znimg_biasfield_upsampled.darr.chunks)
print(znimg_biasfield_upsampled.darr.blocks.shape)


print('biasfield_downsampled')
znimg_biasfield_downsampled =  ZarrNii.from_path(snakemake.params.bf_ds_uri)
print(znimg_biasfield_downsampled.darr.shape)
print(znimg_biasfield_downsampled.darr.chunks)
print(znimg_biasfield_downsampled.darr.blocks.shape)




znimg_biasfield_upsampled.to_ome_zarr(snakemake.params.bf_us_uri,max_layers=0)
print('biasfield_upsampled_input')
znimg_biasfield_upsampled = ZarrNii.from_path(snakemake.params.bf_us_uri)
print(znimg_biasfield_upsampled.darr.shape)
print(znimg_biasfield_upsampled.darr.chunks)
print(znimg_biasfield_upsampled.darr.blocks.shape)




#now multiply biasfield and hires image together,

znimg_hires = ZarrNii.from_path(snakemake.params.spim_uri,channels=[channel_index],level=hires_level,chunks=znimg_biasfield_upsampled.darr.chunks,rechunk=True) #chunk size comes from upsampled array

print('hires_input')
print(znimg_hires.darr.shape)
print(znimg_hires.darr.chunks)
print(znimg_hires.darr.blocks.shape)


znimg_hires.darr = znimg_hires.darr / znimg_biasfield_upsampled.darr

print('n4 corrected before writing')
print(znimg_hires.darr.shape)
print(znimg_hires.darr.chunks)
print(znimg_hires.darr.blocks.shape)



#write the full-res N4-corrected OME-zarr dataset -  this was failing with rechunk problem..
znimg_hires.to_ome_zarr(snakemake.params.spim_n4_uri,max_layers=5)



