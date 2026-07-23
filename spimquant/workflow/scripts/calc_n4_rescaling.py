import nibabel as nib
import numpy as np

uncorr_img = nib.load(snakemake.input.uncorr)
mask_img = nib.load(snakemake.input.mask)
biasfield_img = nib.load(snakemake.input.biasfield)

uncorr = uncorr_img.get_fdata(dtype=np.float32)
mask = mask_img.get_fdata(dtype=np.float32) > 0
biasfield = biasfield_img.get_fdata(dtype=np.float32)

if not mask.any():
    raise ValueError(
        f"Brain mask {snakemake.input.mask} contains no masked voxels; "
        "cannot compute scale and offset."
    )

# Avoid division by zero in bias field
epsilon = np.finfo(np.float32).eps
biasfield_safe = np.maximum(biasfield, epsilon)

# Corrected image = uncorrected / bias field (mimics what N4BiasFieldApply does)
corrected = uncorr / biasfield_safe

# Compute min/max within the masked region
input_min = float(uncorr[mask].min())
input_max = float(uncorr[mask].max())
corrected_min = float(corrected[mask].min())
corrected_max = float(corrected[mask].max())

# Linear mapping: corrected * scale + offset -> [input_min, input_max]
corrected_range = corrected_max - corrected_min
if corrected_range == 0:
    scale = 1.0
    offset = 0.0
else:
    scale = (input_max - input_min) / corrected_range
    offset = input_min - scale * corrected_min

with open(snakemake.output.scale_offset_params, "w") as f:
    f.write(f"scale={scale}\n")
    f.write(f"offset={offset}\n")
