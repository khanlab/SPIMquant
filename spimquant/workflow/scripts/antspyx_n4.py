import ants

antsimg = ants.image_read(snakemake.input.spim)

antsimg_bias = ants.n4_bias_field_correction(antsimg,return_bias_field=True,**snakemake.params.n4_opts)

#write out nifti
ants.image_write(antsimg_bias,snakemake.output.n4_bf_ds)


