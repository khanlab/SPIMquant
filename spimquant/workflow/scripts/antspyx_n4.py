import ants
from zarrnii import ZarrNii
from lib.utils import get_zarr_store, get_channel_index

store = get_zarr_store(snakemake.params.spim_uri)

level=int(snakemake.wildcards.level)
channel_index = get_channel_index(store, snakemake.wildcards.stain)
in_orient = snakemake.config['in_orientation'] #TODO: this is a bit ugly - update the ZarrNii to recognize None  as an orientation, so we can just pass in_orientation
orient_opt = {} if in_orient == None else {'orientation': in_orient}

znimg = ZarrNii.from_ome_zarr(store,level=level,**orient_opt)

print(znimg.darr)
znimg.to_nifti(snakemake.output.spim_ds)

#now perform ants
antsimg = ants.image_read(snakemake.output.spim_ds)

antsimg_bias = ants.n4_bias_field_correction(antsimg,spline_param=(16,16,16),shrink_factor=8,return_bias_field=True)

#write out nifti
ants.image_write(antsimg_bias,snakemake.output.n4_bf_ds)


