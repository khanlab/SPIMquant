import nibabel as nib
nib.load(snakemake.input[0]).to_filename(snakemake.output[0])
