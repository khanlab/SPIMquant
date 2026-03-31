# Output Files Reference

This page describes every output type produced by SPIMquant, what it contains, and where it fits in the processing pipeline.  For a narrative explanation of how each output is generated see [How SPIMquant Works Under the Hood](../workflow_overview.md).

## Directory Structure

All outputs follow [BIDS Derivatives](https://bids-specification.readthedocs.io/en/stable/05-derivatives/01-introduction.html) naming conventions and are written into the `output/spimquant/` root by default.

```
output/spimquant/
├── sub-<label>/
│   ├── micr/                              # Microscopy images and reports
│   │   ├── *_stain-<stain>_level-<N>_SPIM.nii.gz
│   │   ├── *_stain-<stain>_level-<N>_desc-brain_mask.nii.gz
│   │   ├── *_stain-<stain>_level-<N>_desc-N4_SPIM.nii.gz
│   │   ├── *_stain-<stain>_level-<N>_desc-N4_biasfield.nii.gz
│   │   ├── *_stain-<stain>_level-<N>_space-<template>_SPIM.nii.gz
│   │   ├── *_stain-<stain>_level-<N>_space-<template>_regqc.html
│   │   └── *_stain-<stain>_seg-<seg>_from-<template>_level-<N>_patches/
│   ├── parc/                              # Atlas parcellations in subject space
│   │   └── *_seg-<seg>_level-<N>_from-<template>_dseg.nii.gz
│   ├── seg/                               # Segmentation results
│   │   ├── *_stain-<stain>_level-<N>_desc-<desc>_mask.ozx
│   │   ├── *_stain-<stain>_level-<N>_desc-<desc>_fieldfrac.nii.gz
│   │   └── *_stain-<stain>_level-<N>_desc-<desc>_counts.nii.gz
│   ├── tabular/                           # Per-region statistics
│   │   ├── *_seg-<seg>_from-<template>_stain-<stain>_level-<N>_desc-<desc>_regionpropstats.tsv
│   │   ├── *_seg-<seg>_from-<template>_stain-<stain>_level-<N>_desc-<desc>_countstats.tsv
│   │   ├── *_seg-<seg>_from-<template>_desc-<desc>_mergedsegstats.tsv
│   │   └── *_stain-<stain>_desc-<desc>_regionprops.parquet
│   ├── featuremap/                        # Per-region feature maps in template space
│   │   └── *_seg-<seg>_space-<template>_desc-<desc>_<suffix>.nii.gz
│   └── xfm/                              # Registration transforms
│       ├── *_from-subject_to-<template>_xfm.txt
│       ├── *_from-subject_to-<template>_xfm.nii.gz
│       └── *_from-<template>_to-subject_xfm.nii.gz
└── group/                                 # Group-level outputs (no subject subdirectory)
    ├── *_seg-<seg>_from-<template>_desc-<desc>_groupstats.tsv
    └── *_seg-<seg>_from-<template>_desc-<desc>_groupstats.png
```

---

## Participant-Level Outputs

### Downsampled NIfTI Image (`micr/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_SPIM.nii.gz`

A single resolution level extracted from the subject's OME-Zarr multi-scale image and converted to NIfTI format.  The `level` entity indexes into the OME-Zarr pyramid (0 = finest resolution; higher numbers = more downsampled).  `registration_level` in the config controls which level is used for registration; a separate `segmentation_level` is used for segmentation.

**Example from `tests/bids_ds`:**
```
sub-AS36F2/micr/sub-AS36F2_sample-brain_stain-PI_level-3_SPIM.nii.gz
```

---

### Brain Mask (`micr/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_desc-brain_mask.nii.gz`

Binary whole-brain mask (1 = brain, 0 = background) in subject space.  Created by the [masking stage](../workflow_overview.md#stage-3--masking) using Atropos Gaussian mixture modelling combined with a template-space prior.

Used to restrict N4 bias-field correction and registration to voxels within the brain.

---

### Bias-Field Corrected Image (`micr/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_desc-N4_SPIM.nii.gz`

The brain image after N4 (or Gaussian) bias-field correction has been applied.  Intensity non-uniformities across the field of view are substantially reduced, making registration and segmentation more robust.

**Bias field:** `*_stain-<stain>_level-<N>_desc-N4_biasfield.nii.gz` — the multiplicative bias field estimated by N4, saved for QC purposes.

---

### Registered Image in Template Space (`micr/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_space-<template>_SPIM.nii.gz`

The bias-field-corrected subject image warped into the chosen template space (e.g. `space-ABAv3`).  All subjects registered to the same template can be directly compared voxel-by-voxel or displayed as overlays.

This file is the primary input for group-level feature maps.

---

### Registration QC Report (`micr/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_space-<template>_regqc.html`

A self-contained HTML file containing interactive overlays of:
- Subject image in template space (with template contours)
- Template image in subject space (with subject brain-mask contours)

Open this file in any modern web browser to visually verify that the registration is correct before relying on quantitative results.  Common failure modes visible here include flipped axes, large misalignment, or collapsed/inflated brain regions.

---

### Atlas Parcellation in Subject Space (`parc/`)

**Filename pattern:** `*_seg-<seg>_level-<N>_from-<template>_dseg.nii.gz`

Integer-valued parcellation (discrete segmentation) atlas warped from template space into subject space.  Each voxel holds the integer index of the brain region it belongs to; the corresponding region names are in the label TSV file.

This image is used downstream to attribute detected objects and field-fraction values to named brain regions.

---

### Segmentation Mask (`seg/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_desc-<desc>_mask.ozx`

Full-resolution binary segmentation mask stored in compressed OME-Zarr format (`.ozx` extension).  Voxels with value `100` are classified as positive (pathology present); `0` means negative.

!!! note "Scale convention"
    Masks are stored on a 0–100 scale (not 0–1) so that any spatial downsampling of the mask directly yields a *percentage* occupancy (field fraction 0–100 %).

The `<desc>` entity encodes which correction method and segmentation algorithm were used (e.g. `desc-thresholdN4`).

---

### Field-Fraction Map (`seg/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_desc-<desc>_fieldfrac.nii.gz`

A downsampled NIfTI image where each voxel value (0–100) represents the **percentage of high-resolution voxels within that coarser voxel that are positive** according to the segmentation mask.  This converts the binary high-resolution mask into a smooth density map that can be registered to template space.

Field fraction is the primary continuous measure reported per atlas region in the tabular outputs.

---

### Object Count Map (`seg/`)

**Filename pattern:** `*_stain-<stain>_level-<N>_desc-<desc>_counts.nii.gz`

A downsampled NIfTI image where each voxel value is the **number of detected objects** whose centroid falls within that voxel.  Provides a spatial density map of individual pathological objects.

---

### Region-Properties Statistics Table (`tabular/`)

**Filename pattern:** `*_seg-<seg>_from-<template>_stain-<stain>_level-<N>_desc-<desc>_regionpropstats.tsv`

A TSV file with one row per atlas region containing per-object summary statistics:
- Mean, median, and standard deviation of object volume
- Mean, median, and standard deviation of object intensity
- Total object count per region
- Region label index and name (joined from the atlas label TSV)

---

### Count Statistics Table (`tabular/`)

**Filename pattern:** `*_seg-<seg>_from-<template>_stain-<stain>_level-<N>_desc-<desc>_countstats.tsv`

Simpler table with one row per atlas region reporting the raw count of detected objects mapped to that region.

---

### Merged Segmentation Statistics Table (`tabular/`)

**Filename pattern:** `*_seg-<seg>_from-<template>_desc-<desc>_mergedsegstats.tsv`

The primary per-subject per-region statistics file, aggregating metrics across all configured stains.  Column names follow the pattern `<stain>+<metric>` (e.g. `abeta+fieldfrac_mean`, `Iba1+count`).  This is the file consumed by the group-level statistical analysis.

---

### Region-Properties Parquet (`tabular/`)

**Filename pattern:** `*_stain-<stain>_desc-<desc>_regionprops.parquet`

Per-object (not per-region) table in Apache Parquet format containing every detected object and its individual properties (volume, centroid, mean intensity, etc.).  Parquet is used for efficiency — these tables can contain millions of rows for datasets with large numbers of small objects.

---

### Feature Map in Template Space (`featuremap/`)

**Filename pattern:** `*_seg-<seg>_space-<template>_desc-<desc>_<suffix>.nii.gz`

A volumetric NIfTI where each brain region is painted with a scalar value derived from the tabular statistics (e.g. mean field fraction, object count density).  The `<suffix>` encodes the metric name.  These maps are in template space, so any atlas overlay or group comparison is straightforward in a neuroimaging viewer (FSLeyes, ITK-SNAP, etc.).

---

### Registration Transforms (`xfm/`)

**Affine transform:** `*_from-subject_to-<template>_xfm.txt`  
An ITK-compatible affine matrix describing the linear component of the subject-to-template mapping.

**Deformable warp:** `*_from-subject_to-<template>_xfm.nii.gz`  
Dense displacement field for the nonlinear component of the subject-to-template mapping.

**Inverse warp:** `*_from-<template>_to-subject_xfm.nii.gz`  
Inverse displacement field used to warp atlas labels from template space back to each subject.

---

### Image Patches (`micr/`)

**Filename pattern:** `*_stain-<stain>_seg-<seg>_from-<template>_level-<N>_patches/`

A directory of NIfTI files, each a small 3-D crop extracted from a specific atlas region.  Patches are named with the atlas label abbreviation and a patch index (e.g. `CA1_patch-0003.nii.gz`).  Three patch sets are produced:
- `*_patches/` — raw SPIM data
- `*_desc-corrected_patches/` — bias-field-corrected SPIM data
- `*_desc-mask_patches/` — segmentation mask

See [Imaris Crops](../howto/imaris_crops.md) for exporting patches to Imaris format.

---

## Group-Level Outputs

Group-level outputs are produced when running `analysis_level group` and are stored directly under `output/spimquant/group/`.

### Group Statistics Table (`group/`)

**Filename pattern:** `*_seg-<seg>_from-<template>_desc-<desc>_groupstats.tsv`

One row per brain region containing the outcome of statistical tests comparing the two groups defined by `contrast_column` and `contrast_values` in the config.  Columns include:
- `region` — region label name
- `<stain>+<metric>_tstat` — t-statistic
- `<stain>+<metric>_pval` — uncorrected p-value
- `<stain>+<metric>_pval_fdr` — FDR-corrected p-value
- `<stain>+<metric>_effect_size` — Cohen's d or Mann-Whitney effect size

### Group Statistics Heatmap (`group/`)

**Filename pattern:** `*_seg-<seg>_from-<template>_desc-<desc>_groupstats.png`

Annotated heatmap visualising the per-region statistical results.  Rows are brain regions, columns are stain-metric combinations.  Colour encodes effect size; significance markers (stars) indicate FDR-corrected p < 0.05.

---

## Intermediate Files

Many intermediate files are generated during processing and then deleted once the downstream rule has consumed them (marked `temp()` in the Snakemake rules).  These include:
- Pre-Atropos downsampled images
- Gaussian-corrected OME-Zarr (before bias-field-corrected NIfTI is produced)
- Intermediate registration files

If you need to preserve intermediate outputs, run Snakemake with `--notemp`.

---

## Next Steps

- [How SPIMquant Works Under the Hood](../workflow_overview.md) — pipeline narrative
- [Configuration Options](config.md) — controlling output behaviour
- [Group Analysis](../usage/group_analysis.md) — running group-level statistics