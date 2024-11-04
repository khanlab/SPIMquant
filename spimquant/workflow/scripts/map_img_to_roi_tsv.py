
import nibabel as nib
import pandas as pd
import numpy as np

img_nib = nib.load(snakemake.input.img)
dseg_nib = nib.load(snakemake.input.dseg)
dseg_df = pd.read_table(snakemake.input.label_tsv )

img = img_nib.get_fdata()
dseg = dseg_nib.get_fdata()
zooms = dseg_nib.header.get_zooms()

# voxel size in mm^3
voxel_mm3 = np.prod(zooms)

def calc_label_volume(x,dseg,voxel_mm3):
    return np.sum(dseg == x)*voxel_mm3

def calc_avg_fieldfrac(x,img,dseg):
    return np.mean(img[dseg == x])


#calc volume and avg fieldfrac for each label
dseg_df['volume'] = dseg_df['index'].apply(calc_label_volume,args=(img,voxel_mm3))
dseg_df['avg_fieldfrac'] = dseg_df['index'].apply(calc_avg_fieldfrac,args=(img,dseg))


dseg_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

