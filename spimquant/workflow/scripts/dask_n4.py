from zarrnii import ZarrNii
from lib.utils import get_zarr_store, get_channel_index


store = get_zarr_store(snakemake.params.spim_uri)

channel_index = get_channel_index(store, snakemake.wildcards.stain)

in_orient = snakemake.config['in_orientation'] #TODO: this is a bit ugly - update the ZarrNii to recognize None  as an orientation, so we can just pass in_orientation
orient_opt = {} if in_orient == None else {'orientation': in_orient}

hires_level=int(snakemake.wildcards.level)
ds_level=int(snakemake.wildcards.dslevel)

#this function optionally uses coiled


hires_shape = ZarrNii.from_ome_zarr(store,channels=[channel_index],level=hires_level,**orient_opt).darr.shape

#ok, now we have the bias field.. let's resample it to the level
# where we perform thresholding

# Compute step 1: first, we write the nifti n4 bias field to ome zarr (this is so we can use distributed computing)
ZarrNii.from_nifti(snakemake.input.n4_bf_ds,chunks=(1,10,10,10)).to_ome_zarr(snakemake.params.bf_ds_uri)

if snakemake.config['use_coiled']:
    from coiled import Cluster
    cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=10)
    client = cluster.get_client()


# Compute step 2: then we upsample the bias field, and save to file
znimg_biasfield_upsampled = ZarrNii.from_ome_zarr(snakemake.params.bf_ds_uri,**orient_opt).upsample(to_shape=hires_shape)
znimg_biasfield_upsampled.to_ome_zarr(snakemake.params.bf_us_uri,max_layer=0)


# Compute step 3: now multiply biasfield and hires image together.  we read the biasfield 
# in lazily again (this is so the subsequent step is a different dask computation than the previous step), 
# and write the full-res N4-corrected OME-zarr dataset
znimg_biasfield_upsampled = ZarrNii.from_ome_zarr(snakemake.params.bf_us_uri,**orient_opt)
znimg_hires = ZarrNii.from_ome_zarr(store,channels=[channel_index],level=hires_level,chunks=znimg_biasfield_upsampled.darr.chunks,rechunk=True,**orient_opt) #chunk size comes from upsampled array
znimg_hires.darr = znimg_hires.darr / znimg_biasfield_upsampled.darr
znimg_hires.to_ome_zarr(snakemake.params.spim_n4_uri,max_layer=5)


