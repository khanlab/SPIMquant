# Reusing Outputs from a Previous Run

Brain masking and template registration are the same regardless of how you
segment and quantify afterwards. With `--reuse_derivatives` you can compute
them once and then point later runs at that output dataset, instead of
recomputing them for every new segmentation configuration.

## Example

```bash
# 1. Register only: produces brain masks and transforms
spimquant /path/to/bids /path/to/derivatives/spimquant_reg participant \
    --no_segmentation --no_vessels --cores all

# 2. Segment in a new output dataset, reusing the masks and transforms
spimquant /path/to/bids /path/to/derivatives/spimquant_gmm participant \
    --reuse_derivatives /path/to/derivatives/spimquant_reg \
    --seg_method gmm+n3k1 --cores all
```

At startup SPIMquant reports how many subjects it found each output for, e.g.

```
Reusing brain masks from /path/to/derivatives/spimquant_reg for 12 of 12 subject(s); the rest will be computed
```

## What is reused

Select with `--reuse_outputs` (default: both):

| Value        | Files                                                                      |
|--------------|----------------------------------------------------------------------------|
| `brainmask`  | `micr/*_stain-<reg stain>_level-<registration_level>_desc-brain_mask.nii.gz` |
| `transforms` | the composite warps `xfm/*_from-subject_to-<template>_xfm.nii.gz` / `*_from-<template>_to-subject_xfm.nii.gz` |

Reused files are copied into the new output directory, so it stays a complete,
self-contained dataset. To save disk space (the composite warps are ~1.5 GB per
subject), pass `--reuse_symlink` to symlink them instead; the new dataset then
depends on the prior one staying in place.

## Matching and fallback

Files are looked up at the exact BIDS paths SPIMquant itself writes for the
current settings, so they are only reused when the path-defining settings match:
the subject/session/sample/acquisition entities, `--template`,
`--registration_level` and the registration stain. Anything not found is
computed as usual, per subject and per file, including subjects missing from
the prior dataset.

The registration QC report (`xfm/*_regqc.html`) is not made for subjects whose
transforms are reused, since it would only re-check the prior run; see the
prior dataset's report instead.

Settings that change the result without changing the filenames, such as
`--orientation`, `--template_crop` or `--sloppy`, are **not** checked. Use the
same values as the run that produced the dataset you are reusing.
