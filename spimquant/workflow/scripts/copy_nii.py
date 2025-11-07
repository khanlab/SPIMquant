import nibabel as nib
img = nib.load(snakemake.input[0]).to_filename(snakemake.output[0])
