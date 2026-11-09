"""Template-space fibre metric maps from warped moments.

The moments (a 4D NIfTI, one volume per moment) are warped to the template
with linear interpolation, and the metrics are computed afterwards.  Since
every metric is a ratio of moments (e.g. fibre voxels / tissue voxels), and
the dispersion and anisotropy metrics depend only on rotation-invariant
quantities (tensor eigenvalues, SH power per order), this is valid without
reorienting the tensors or ODFs, and avoids interpolating across voxels
with and without fibres.  Direction-dependent outputs (principal direction,
ODF) are not produced here: they would need reorientation by the warp.

Only the non-SH moments are warped (warping all SH sums onto a template
grid takes tens of GB), so the ODF GFA is not available in template space;
it is in the subject-space maps and the atlas-region tables.
"""

import numpy as np


def main():
    import nibabel as nib
    from fibre_orientation import N_FIBRE, N_TISSUE, SH_START, moments_to_metrics

    img = nib.load(snakemake.input.moments_nii)
    data = np.asarray(img.dataobj, dtype=np.float64)  # x, y, z, moment
    moments = np.moveaxis(data, -1, 0)[:SH_START]
    # linear interpolation cannot make counts negative, but guard round-off
    moments[:2] = np.clip(moments[:2], 0, None)
    # a single l=0 SH channel (its value is irrelevant to the metrics kept)
    moments = np.concatenate([moments, moments[N_FIBRE][None]])
    metrics = moments_to_metrics(moments)
    # Near the mask edge the interpolated counts become tiny and their ratios
    # are noise: keep voxels with at least min_tissue_fraction of a full
    # super-voxel of tissue (a full one taken as the 99th percentile).
    tissue = moments[N_TISSUE]
    full = np.percentile(tissue[tissue > 0], 99) if np.any(tissue > 0) else 0.0
    keep = tissue >= float(snakemake.params.min_tissue_fraction) * full
    for name in snakemake.params.metrics:
        data = np.where(keep, metrics[name], 0).astype(np.float32)
        out = nib.Nifti1Image(data, img.affine)
        out.header.set_xyzt_units("mm")
        nib.save(out, getattr(snakemake.output, name))


if __name__ == "__main__":
    main()
