# Multiple scans of a sample

A sample is sometimes imaged more than once, e.g. a whole-brain scan with a 1x
objective and a higher-resolution scan with a 12x objective that covers only part
of the brain. SPIMquant can process both: one scan per sample, the **reference
scan**, is registered to the template as usual, and each **additional scan** is
registered to the reference scan of its sample. The additional scans get atlas
labels, template-space outputs and quantification through the reference scan's
registration to the template, so their field of view does not need to hold
enough of the brain to register to the template on its own.

## Selecting the scans

The reference scans come from the `spim` input and the additional scans from the
`spim_aux` input, which is only used with `--aux_spim`. Both select every SPIM
scan by default, so use `--filter-spim` and `--filter-spim_aux` to split them,
typically by the `acq-` entity:

```bash
spimquant /path/to/bids /path/to/output participant \
  --aux_spim \
  --filter-spim acquisition=1x \
  --filter-spim_aux acquisition=12x
```

for a dataset such as:

```
sub-01/micr/sub-01_sample-brain_acq-1x_SPIM.ome.zarr
sub-01/micr/sub-01_sample-brain_acq-12x_SPIM.ome.zarr
sub-02/micr/sub-02_sample-brain_acq-1x_SPIM.ome.zarr
sub-02/micr/sub-02_sample-brain_acq-12x_SPIM.ome.zarr
```

Each additional scan is paired with the reference scan that has the same
`subject`, `session` and `sample` entities. An additional scan without a
reference scan is skipped with a warning, and it is an error for a sample to have
more than one reference scan or for a scan to be selected as both.

The reference and additional scans must have the same BIDS entities (e.g. both an
`acq-`), and the same stains (channels).

## Registration to the reference scan

The additional scan is registered with greedy, using all stains together
(multichannel; `--aux_reg_stains` to choose them), at `--registration_level` of
each scan. The additional scan is the fixed image, so the similarity metric is
computed over its field of view.

`--aux_reg_method` sets the transformation:

- `rigid`: for scans of the same mount with the same calibration
- `affine` (default): rigid, then affine, which also takes up scaling differences
  between objectives (e.g. axial scaling from a refractive index mismatch)
- `deform`: affine, then deformable, e.g. for a sample re-mounted between scans

The physical coordinates of two scans need not relate to each other (each scan's
origin is typically wherever its stage started), so by default
(`--aux_reg_init match`) the registration starts from the translation that best
places the additional scan within the reference scan, found by normalized
cross-correlation over every position. This assumes the two scans share their
axes, as they do in the same mount. If they do share a stage coordinate system,
`--aux_reg_init physical` starts from the physical coordinates instead, and
`--aux_reg_init centers` aligns the image centers.

The metric, iterations and smoothing are in the `aux_reg` section of the config.

## Outputs

For an additional scan (e.g. `acq-12x`), outputs have the same names as for any
scan, with its own entities. In addition:

- `xfm/..._acq-12x_from-ref_to-subject_desc-affine_type-ras_xfm.txt`: the greedy
  affine from the additional scan to its reference scan
- `xfm/..._acq-12x_space-subject_from-ref_stain-<stain>_SPIM.nii.gz`: the
  reference scan resampled into the additional scan, to check the registration
  (overlay it on the additional scan's `level-<registration_level>` image)
- `xfm/..._acq-12x_from-subject_to-<template>_desc-viaref_xfm.nii.gz` (and its
  inverse): the composite warps to and from the template, through the reference
  scan

The brain mask of an additional scan uses the reference scan's brain mask,
resampled into it, as its spatial prior (in place of the template mask), so it is
only as good as the reference scan's.

Only reference scans are used for template registration QC reports, `tune_reg`,
cohort templates, MRI registration and the `group` analysis level (where the
additional scans would count their participants twice); the additional scans'
per-scan tables are under each scan's own entities.

!!! note
    Atlas regions only partly covered by an additional scan's field of view have
    their volume (and so density) computed over the covered part only.
