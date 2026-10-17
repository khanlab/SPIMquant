import nibabel as nib
import numpy as np
from zarrnii import ZarrNii

target_level = int(snakemake.wildcards.level)
hires_level = int(snakemake.params.hires_level)

downsampling_level = target_level - hires_level
if downsampling_level < 0:
    raise ValueError("Target level for fieldfrac is smaller than the input level!")


# this will downsample automatically based on the level
znimg_density_ds = ZarrNii.from_file(
    snakemake.input.mask,
    level=downsampling_level,
    downsample_near_isotropic=True,
)

print(znimg_density_ds.darr)
fieldfrac_nii = znimg_density_ds.to_nifti()
ref_nii = nib.load(snakemake.input.ref)

# Put the field fraction on the grid of the SPIM image at the same level, which
# the atlas parcellation is also on. Both are block averages starting from voxel 0
# of the same full-resolution image, so voxel i covers the same tissue in each,
# but they disagree on the rest of the grid: on an odd dimension the mask pyramid
# rounds the size down where the SPIM pyramid rounds it up, and the mask pyramid
# levels lack the half-voxel translation of a block average. So keep the data,
# indexed from voxel 0, and take the shape and affine from the reference. A
# trailing partial block missing from the field fraction is filled with zeros.
fieldfrac = np.asanyarray(fieldfrac_nii.dataobj).squeeze()
ref_shape = ref_nii.shape[:3]

if not np.allclose(fieldfrac_nii.affine[:3, :3], ref_nii.affine[:3, :3], rtol=1e-3):
    raise ValueError(
        "Field fraction voxel axes do not match the reference image:\n"
        f"{fieldfrac_nii.affine}\nvs\n{ref_nii.affine}"
    )
if fieldfrac.ndim != 3 or any(
    abs(n - n_ref) > 1 for n, n_ref in zip(fieldfrac.shape, ref_shape)
):
    raise ValueError(
        f"Field fraction shape {fieldfrac.shape} is not within one voxel of the "
        f"reference shape {ref_shape}"
    )

fieldfrac_on_ref = np.zeros(ref_shape, dtype=fieldfrac.dtype)
overlap = tuple(slice(0, min(n, n_ref)) for n, n_ref in zip(fieldfrac.shape, ref_shape))
fieldfrac_on_ref[overlap] = fieldfrac[overlap]

out_nii = nib.Nifti1Image(fieldfrac_on_ref, ref_nii.affine)
out_nii.header.set_xyzt_units("mm")
out_nii.to_filename(snakemake.output.fieldfrac_nii)
