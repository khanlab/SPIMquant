import nibabel as nib
import numpy as np
import pandas as pd

dseg_df = pd.read_table(snakemake.input.label_tsv)


img_nib = nib.load(snakemake.input.dseg)
img = img_nib.get_fdata()
zooms = img_nib.header.get_zooms()

# voxel size in mm^3
voxel_mm3 = np.prod(zooms)


def calc_label_volume(x, img, voxel_mm3):
    return np.sum(img == x) * voxel_mm3


dseg_df["volume"] = dseg_df["index"].apply(calc_label_volume, args=(img, voxel_mm3))

dseg_df.to_csv(snakemake.output.volumes_tsv, sep="\t", index=False)
