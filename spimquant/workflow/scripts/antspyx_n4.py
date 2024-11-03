import ants
from zarrnii import ZarrNii


level=int(snakemake.wildcards.level)
level_z=level-1 #z downsampling one less (since already lower-res

ZarrNii.from_path(snakemake.params.spim_uri,level=level).downsample(along_z=level_z).to_nifti(snakemake.output.spim_ds)

#now perform ants
antsimg = ants.image_read(snakemake.output.spim_ds)

antsimg_bias = ants.n4_bias_field_correction(antsimg,spline_param=(16,16,16),shrink_factor=8,return_bias_field=True)

#write out nifti
ants.image_write(antsimg_bias,snakemake.output.n4_bf_ds)


