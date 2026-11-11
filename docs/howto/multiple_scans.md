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
  reference scan resampled into the additional scan
- `xfm/..._acq-12x_from-ref_to-subject_regqc.html`: a QC report of the
  registration, with the reference scan resampled into the additional scan with the
  initial translation (and its match score) and with the registration, for each of
  `--aux_reg_stains`, and the atlas labels the additional scan gets through it
- `xfm/..._acq-12x_from-subject_to-<template>_desc-viaref_xfm.nii.gz` (and its
  inverse): the composite warps to and from the template, through the reference
  scan
- a `coverage` column in `tabular/..._acq-12x_..._mergedsegstats.tsv`: the
  fraction (0–1) of each atlas region, as labelled in the reference scan, that the
  additional scan's field of view (its nonzero voxels) covers

The brain mask of an additional scan uses the reference scan's brain mask,
resampled into it, as its spatial prior (in place of the template mask), so it is
only as good as the reference scan's.

Only reference scans are used for template registration QC reports, `tune_reg`,
cohort templates and MRI registration.

!!! note
    The volume, and so the density, of an atlas region in an additional scan is over
    the part its field of view covers. Field fraction and density are comparable
    between scans when the covered part is representative of the region; volumes
    and counts are not, so use `coverage` to judge them (the covered volume is
    about `coverage` times the region's volume in the reference scan).

## Group analysis

At the `group` analysis level, the reference scans give the usual outputs, and
the additional scans give their own for each acquisition (the entities that tell
the scans of a sample apart, e.g. `acq-`), named with it, so that no participant
counts twice in the same statistics:

```
group/<label>/seg-roi22_from-ABAv3_desc-otsu+k3i2_allsubjects.tsv
group/<label>/acq-12x_seg-roi22_from-ABAv3_desc-otsu+k3i2_allsubjects.tsv
group/<label>/acq-12x_seg-roi22_from-ABAv3_..._contrast-<contrast>_..._groupstats.tsv
```

The additional scans have the merged tables (`allsubjects.tsv`, with the coverage
column), and with `--group-stats-pairwise`, the group statistics, their heatmaps and
template-space maps. In their statistics, a region's metrics are only used for
scans that cover at least `--aux_min_coverage` (default 0.5) of it, so the number
of scans in the statistics can differ between regions. The template-space count maps and the
regionprops tables of the group level are for the reference scans only.
