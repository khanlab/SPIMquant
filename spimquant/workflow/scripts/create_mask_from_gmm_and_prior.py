import nibabel as nib
import numpy as np

tissue_nib = nib.load(snakemake.input.tissue_dseg)
template_mask_nib = nib.load(snakemake.input.template_mask)

tissue_vol = tissue_nib.get_fdata()
template_mask_vol = template_mask_nib.get_fdata()

# now, see what labels are in the template_mask:
fg_tissue = tissue_vol * (template_mask_vol > 0)
bg_tissue = tissue_vol * (template_mask_vol == 0)

out_mask = np.zeros(template_mask_vol.shape)

for i in range(1, snakemake.params.k + 1):

    # if more voxels in foreground than in background, we assign it to the mask
    nvox_fg = (fg_tissue == i).sum()
    nvox_bg = (bg_tissue == i).sum()
    print(f"i={i}, {nvox_fg} to {nvox_bg} (fg to bg)")
    if nvox_fg > nvox_bg:

        out_mask[tissue_vol == i] = 1

out_nib = nib.Nifti1Image(out_mask, affine=tissue_nib.affine, header=tissue_nib.header)
out_nib.to_filename(snakemake.output.mask)
