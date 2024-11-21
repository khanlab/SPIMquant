import nibabel as nib
import numpy as np
from skimage.filters import threshold_multiotsu

vol = nib.load(snakemake.input.masked).get_fdata()


# Mask the zero-valued voxels
nonzero_mask = vol > 0
nonzero_values = vol[nonzero_mask]

multi_thresholds = threshold_multiotsu(nonzero_values,classes=snakemake.params.otsu_n_classes)

np.save(snakemake.output.otsu_thresholds,multi_thresholds)
