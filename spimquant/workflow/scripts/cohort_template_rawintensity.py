"""Voxelwise mean of subject images that are already in cohort template space."""

import nibabel as nib
import numpy as np

total = 0
for path in snakemake.input.images:
    img = nib.load(path)
    total = total + np.asarray(img.dataobj, dtype=np.float64)

out_img = nib.Nifti1Image(
    (total / len(snakemake.input.images)).astype(np.float32), img.affine, img.header
)
out_img.set_data_dtype(np.float32)
out_img.to_filename(snakemake.output[0])
