import nibabel as nib
import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input.tsv, sep="\t")
dseg_nib = nib.load(snakemake.input.dseg)
dseg_vol = dseg_nib.get_fdata()

out_vol = np.zeros(dseg_vol.shape)

for i in range(len(df)):
    label = df.loc[i, snakemake.params.label_column]
    feature = df.loc[i, snakemake.params.feature_column]

    out_vol[dseg_vol == label] = feature

out_nib = nib.Nifti1Image(out_vol, affine=dseg_nib.affine)
out_nib.to_filename(snakemake.output.nii)
