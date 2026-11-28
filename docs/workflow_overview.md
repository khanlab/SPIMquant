# How it works

This page follows one brain through the `participant` level, then the
`group` level. The [workflow graphs](reference/workflow.md) show every rule
and its dependencies, and the [output reference](reference/outputs.md) every
file; both are generated from the workflow.

Two resolutions are used throughout. Registration, masking and all per-region
maps work on a coarse level of the image pyramid (`--registration-level`,
default 5: about 50 µm in-plane for a scan at 1.6 µm). Segmentation works at full
resolution (`--segmentation-level`, default 0), and its results are brought
back to the coarse grid as fractions and counts.

## 1. Brain mask

The registration stain (the first of `--stains-for-reg` in the scan) is read
at the registration level. A Gaussian mixture (ANTs Atropos) separates
tissue from background, and the template's brain mask, brought in by a quick
affine registration, removes what is not brain. The mask (`micr/*_desc-brain_mask`)
limits bias correction, registration and segmentation to the brain.

## 2. Registration

The registration image is N4 bias-corrected and masked, then registered to
the template: affine, then deformable, with the method set by
`--registration_method` (default `greedy`; also ANTs and FireANTs presets,
see [Registration methods](howto/registration_methods.md)). The affine and
warp are composed into one warp in each direction (`xfm/`). Every later step
that moves data between subject and template uses these two files.

For cohorts whose contrast differs from the template's, the subjects can be
registered to a [cohort template](howto/cohort_template.md) built from
themselves, whose single registration to the atlas is then reused.

The registration QC report (`xfm/*_regqc.html`) overlays the atlas on the
warped subject and the subject on the template.

## 3. Atlas labels in the subject

The template's atlases (`--atlas_segs`) are warped into the subject
(`parc/*_dseg.nii.gz`). Objects and voxels are assigned to regions in subject
space, so no full-resolution data is ever resampled for quantification.

## 4. Intensity correction

Each segmented stain gets its own bias field, fitted at
`--correction-level` (default 6) with N4 (`--correction-method`, default
`n4`; also `gaussian`, or `none`/`log` for methods that respond to local
contrast) and divided out of the full-resolution image.

## 5. Segmentation

Each stain is segmented with its mask method and split into objects with its
object method, both set per stain (`--seg_method`, `--object_method`; see
[Segmentation](howto/segmentation.md)):

- histogram methods (`gmm+...`, the default; `otsu+...`; `th<value>`) pick a
  global threshold from the corrected image of the whole brain;
- `dog+...` detects small puncta with a difference of Gaussians;
- `lantern` runs a deep-learning ensemble for amyloid-beta plaques on a
  GPU (the default for Abeta when one is available);
- vessels (CD31, lectin) are segmented with vesselFM, and their signed
  distance transform is kept so that every object gets its distance to the
  nearest vessel.

Masks are stored at full resolution as OME-Zarr with values 0 and 100, so
that downsampling them gives the percent of each coarse voxel that is
covered: the **field fraction** (`seg/*_fieldfrac.nii.gz`).

Objects are connected components of the mask (or a seeded watershed, which
separates touching cells). For each object SPIMquant records its centroid in
subject and template space, size, intensity and shape. It then finds, for
every pair of stains, the distance from each object to the nearest object of
the other stain, and which objects overlap (colocalization).

## 6. Per-region tables and maps

Field fraction, object counts, densities and object properties are summed or
averaged within each atlas region:

- `tabular/*_segstats.tsv`: one stain and one choice of methods;
- `tabular/*_mergedsegstats.tsv`: all stains (primary methods),
  colocalization and fibre metrics in one row per region. This is the table
  the group level uses;
- `tabular/*_regionprops.parquet`: every object of every stain, with its
  distances, for analyses beyond region means;
- `featuremap/`: each region's value painted into the template, one image per
  stain and metric, ready to view over the atlas or average across brains.

Field fraction and count maps are also warped to the template voxel by voxel
(`seg/*_space-{template}_*`).

QC figures for each stain go to `qc/`: segmentation overlays on slices and
zoomed regions, intensity histograms with the chosen threshold, object size
distributions and depth profiles.

## 7. Group level

`spimquant BIDS OUT group` reads the participant outputs in `OUT`:

- stacks every subject's `mergedsegstats.tsv` and `regionpropstats.tsv`,
  joined with `participants.tsv` (`group/<label>/*_allsubjects.tsv`), and the
  object and colocalization tables of all subjects;
- with `--group-stats-model`, fits an ordinary least-squares model per region
  and metric and, for each pair of levels of the `--group-stats-pairwise`
  factors, writes t-statistics, p-values (uncorrected) and Cohen's d as a
  table, a heatmap and template-space maps. See
  [Group analysis](usage/group_analysis.md).

## Optional analyses

| Option | Adds |
|---|---|
| `--fibres` | Fibre orientation (multiscale Frangi filter, orientation distributions, dispersion and anisotropy maps, per-region fibre metrics) |
| `--register_to_mri` | Registration of the lightsheet scan to the same brain's MRI, and of both to the template ([MRI registration](howto/mri_registration.md)) |
| `--aux_spim` | Processing of extra scans of the same brain, registered to the reference scan ([Multiple scans](howto/multiple_scans.md)) |
| `--template-spim-res`, `--bbox_config`, `--save_patches`, `--crop_labels` | Full-resolution data in template space, crops and patches ([Template-space outputs](howto/template_space_outputs.md)) |
| `--reuse-derivatives` | Reuse of masks and transforms from an earlier run ([Reusing outputs](howto/reuse_derivatives.md)) |
