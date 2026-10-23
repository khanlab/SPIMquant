"""Prepare subject images for FireANTs cohort template building.

Each subject image (already in the original template space) is clipped to
robust percentiles of its nonzero voxels and rescaled to [0, 1]. Writes:

- the normalized images, one per subject and stain
- an image list with one comma-separated row per subject (one entry per stain,
  in the stain order of the template channels)
- the FireANTs build_template config
"""

import os

import nibabel as nib
import numpy as np
import yaml

rows = snakemake.params.rows
stains = snakemake.params.stains
p_lo, p_hi = snakemake.params.clip_percentiles
out_dir = snakemake.output.images

os.makedirs(out_dir, exist_ok=True)

image_list_rows = []
for row in rows:
    out_paths = []
    for stain, in_path in zip(stains, row["images"]):
        img = nib.load(in_path)
        data = np.asarray(img.dataobj, dtype=np.float32)
        nonzero = data[data > 0]
        if nonzero.size == 0:
            raise ValueError(f"{in_path} has no nonzero voxels")
        lo, hi = np.percentile(nonzero, [p_lo, p_hi])
        data = (np.clip(data, lo, hi) - lo) / max(hi - lo, np.finfo(np.float32).eps)

        out_path = os.path.abspath(
            os.path.join(out_dir, f"{row['name']}_stain-{stain}.nii.gz")
        )
        out_img = nib.Nifti1Image(data, img.affine, img.header)
        out_img.set_data_dtype(np.float32)
        out_img.to_filename(out_path)
        out_paths.append(out_path)
        print(f"{in_path}: clipped to [{lo:.4g}, {hi:.4g}] -> {out_path}")
    image_list_rows.append(",".join(out_paths))

with open(snakemake.output.image_list, "w") as f:
    f.write("\n".join(image_list_rows) + "\n")

os.makedirs(os.path.dirname(snakemake.output.fireants_config), exist_ok=True)
with open(snakemake.output.fireants_config, "w") as f:
    yaml.safe_dump(dict(snakemake.params.fireants), f, sort_keys=False)
