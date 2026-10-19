"""
Create a brain mask by combining GMM tissue classes with a template-derived prior.

Rather than assigning each GMM class wholesale to brain or background (which
fails when a class mixes brain and background, e.g. agar-embedded samples
where the background is not the darkest class), each voxel gets a posterior
probability of brain from:

  - a spatial prior: the affine-warped template mask, blurred to allow for
    registration error, and
  - an intensity likelihood ratio per GMM class: how often the class occurs in
    the confident brain core (eroded prior) vs confident background (outside
    the dilated prior).

posterior odds = prior odds * likelihood ratio

A class that is pure background (air or agar) overrides the prior, a class that
is pure brain fills in where the prior is too tight, and a mixed class defers
to the prior. The result is thresholded at 0.5, restricted to voxels that were
classified (i.e. non-zero in the image), reduced to the largest connected
component and hole-filled. Finally, the mask is grown (geodesic dilation, up to
grow_mm) into voxels whose class strongly favours brain (LR >= grow_min_lr), to
recover tissue the affine prior missed without leaking into mixed classes.
"""

import nibabel as nib
import numpy as np
from scipy import ndimage

prior_sigma_mm = snakemake.params.prior_sigma_mm
margin_mm = snakemake.params.margin_mm
grow_mm = snakemake.params.grow_mm
grow_min_lr = snakemake.params.grow_min_lr
eps = 1e-3

tissue_nib = nib.load(snakemake.input.tissue_dseg)
template_mask_nib = nib.load(snakemake.input.template_mask)

tissue_vol = np.asarray(tissue_nib.dataobj).astype(int)
prior_mask = np.asarray(template_mask_nib.dataobj) > 0
zooms = np.array(tissue_nib.header.get_zooms()[:3], dtype=float)

classified = tissue_vol > 0

# confident brain core and confident background, for estimating class likelihoods
margin_vox = np.maximum(np.round(margin_mm / zooms).astype(int), 1)
struct = np.ones(2 * margin_vox + 1, dtype=bool)
core = ndimage.binary_erosion(prior_mask, structure=struct) & classified
outside = ~ndimage.binary_dilation(prior_mask, structure=struct) & classified

n_core = max(core.sum(), 1)
n_outside = max(outside.sum(), 1)
print(f"core voxels: {n_core}, outside voxels: {n_outside}")

# per-class log likelihood ratio, with add-one smoothing so empty classes stay finite
tissue_labels = np.unique(tissue_vol[classified])
log_lr = np.zeros(tissue_vol.max() + 1)
for i in tissue_labels:
    p_brain = ((tissue_vol[core] == i).sum() + 1) / (n_core + len(tissue_labels))
    p_bg = ((tissue_vol[outside] == i).sum() + 1) / (n_outside + len(tissue_labels))
    log_lr[i] = np.log(p_brain / p_bg)
    print(f"i={i}, P(class|brain)={p_brain:.4f}, P(class|bg)={p_bg:.4f}, LR={np.exp(log_lr[i]):.3g}")

# spatial prior: blurred template mask
spatial_prior = ndimage.gaussian_filter(
    prior_mask.astype(np.float32), sigma=prior_sigma_mm / zooms
)
spatial_prior = np.clip(spatial_prior, eps, 1 - eps)

log_odds = np.log(spatial_prior / (1 - spatial_prior)) + log_lr[tissue_vol]
out_mask = (log_odds > 0) & classified

# keep largest connected component
labelled, n_comp = ndimage.label(out_mask)
if n_comp > 1:
    sizes = ndimage.sum_labels(out_mask, labelled, index=np.arange(1, n_comp + 1))
    out_mask = labelled == (np.argmax(sizes) + 1)

# grow into strongly brain-like voxels connected to the mask
grow_iters = int(np.round(grow_mm / zooms.min()))
if grow_iters > 0:
    candidates = out_mask | (log_lr[tissue_vol] >= np.log(grow_min_lr)) & classified
    out_mask = ndimage.binary_dilation(
        out_mask, iterations=grow_iters, mask=candidates
    )

out_mask = ndimage.binary_fill_holes(out_mask) & classified

print(f"mask fraction: {out_mask.mean():.3f} (prior fraction: {prior_mask.mean():.3f})")

out_nib = nib.Nifti1Image(
    out_mask.astype(np.uint8), affine=tissue_nib.affine, header=tissue_nib.header
)
out_nib.set_data_dtype(np.uint8)
out_nib.to_filename(snakemake.output.mask)
