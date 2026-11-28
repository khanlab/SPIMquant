# Reusing Outputs from a Previous Run

Brain masking and template registration are the same regardless of how you
segment and quantify afterwards. With `--reuse_derivatives` you can compute
them once and then point later runs at that output dataset, instead of
recomputing them for every new segmentation configuration.

Segmentation masks and vessel masks can be reused the same way (opt-in), for
runs that only change what is computed from them: object methods, atlases,
tables, group statistics. These are the GPU (LANTERN, vesselFM) and
high-memory (signed distance transform) steps of a run.

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

To also reuse the segmentations and vessels of an earlier full run, e.g. to
quantify with a different object method or atlas:

```bash
spimquant /path/to/bids /path/to/derivatives/spimquant_watershed participant \
    --reuse_derivatives /path/to/derivatives/spimquant_gmm \
    --reuse_outputs brainmask transforms segmentation vessels \
    --seg_method gmm+n3k1 --object_method watershed --cores all
```

## What is reused

Select with `--reuse_outputs` (default: `brainmask transforms`):

| Value          | Files                                                                      |
|----------------|----------------------------------------------------------------------------|
| `brainmask`    | `micr/*_stain-<reg stain>_level-<registration_level>_desc-brain_mask.nii.gz` |
| `transforms`   | the composite warps of `--registration_method`: `xfm/*_from-subject_to-<template>_desc-<method>_xfm.nii.gz` / `*_from-<template>_to-subject_desc-<method>_xfm.nii.gz` |
| `segmentation` | for every `--seg_method` method of every stain (primary and comparison): `seg/*_stain-<stain>_level-<segmentation_level>_desc-<method>_mask.ozx`, with the `_thresholds.png` of gmm, otsu and dog methods and the `_poolmask.nii.gz` of dog; the LANTERN vote map `seg/*_level-neariso<res>_desc-lantern_probseg.ozx`; the signed distance transforms `seg/*_desc-<method>_dist.ozx` of stains in the `sdt_stains` config key |
| `vessels`      | `vessels/*_stain-<vessel stain>_level-<segmentation_level>_desc-vesselfm_mask.ozx` and its signed distance transform `_dist.ozx`; if present, the skeleton `desc-vesselfm+skeleton_mask.ozx` and `_graph`, `_nodes`, `_edges.parquet` (`--vessel_graph`) |

`segmentation` and `vessels` are not on by default: with them on, a run that
changes a segmentation setting the file names do not show (see below) would
not segment again. Everything made from the masks (objects, distances,
colocalization, field fractions, tables, QC) is computed in the new run.

Reused files are copied into the new output directory, so it stays a complete,
self-contained dataset. To save disk space (the composite warps are ~1.5 GB per
subject, and level-0 masks and distance transforms are larger), pass
`--reuse_symlink` to symlink them instead; the new dataset then depends on the
prior one staying in place. The zarr format must match: an `.ome.zarr` mask is
not reused by a run with `--zarr-store-ext ozx`, or the other way round.

## Matching and fallback

Files are looked up at the exact BIDS paths SPIMquant itself writes for the
current settings, so they are only reused when the path-defining settings match:
the subject/session/sample/acquisition entities, `--template`,
`--registration_level`, the registration stain and `--registration_method`.
Datasets made before registration methods were added have warps without a `desc-`;
these were computed with greedy, so they are reused for `--registration_method
greedy` only. They may also have been registered to a packaged SPIM template, which
is not accurately aligned to the atlases, so SPIMquant warns about them: register
again if you need accurate labels. Anything not found is
computed as usual, per subject and per file, including subjects missing from
the prior dataset.

With `--reg-cohort-template`, the `desc-viacohort` warps are reused, and subjects
without them are registered to the cohort template of the prior dataset
(`cohort_template/tpl-<name>/`) rather than a new template being built; see
[Cohort templates](cohort_template.md#in-one-run-reg-cohort-template).

With `--via-template`, the `desc-via<name>` warps are reused, or else the
subject's warps to the via template, with only the rest registered.

The registration QC report (`xfm/*_regqc.html`) is not made for subjects whose
transforms are reused, since it would only re-check the prior run; see the
prior dataset's report instead.

For brain masks and transforms, settings that change the result without
changing the filenames, such as `--orientation`, `--template_crop` or
`--sloppy`, are **not** checked. Use the same values as the run that produced
the dataset you are reusing.

## Segmentation and vessel settings

The file name of a mask gives its stain, level and method, but not every
setting it depends on. Before reusing these outputs, SPIMquant compares the
settings below with the prior dataset's `config/snakebids.yml` (the resolved
configuration of the latest run that wrote to it). If any differs, or the prior
configuration does not have it, that group is computed instead of reused, and
the startup message names the setting:

```
Not reusing th/otsu/gmm segmentation masks from /path/to/derivatives/spimquant_gmm, settings differ: correction_method: prior 'gaussian', now 'n4'. They will be computed.
```

| Group | Settings compared |
|---|---|
| th, otsu and gmm masks | `--correction_method`, `--correction_level`, `--correction_n4_iters`, `--correction_n4_spline_spacing`, `--seg_hist_bin_width`, `--seg_hist_percentile_range` |
| dog masks | the correction settings above, and the `dog_pool_*` and `dog_stat_*` config keys |
| LANTERN vote maps | `--lantern_config`, `plaque_level`, `plaque_iso_res`, `plaque_n_folds`, `plaque_tile`, `plaque_stride`, `--merge_mode` (`plaque_merge_mode`), the model URL (`models.lantern`) |
| LANTERN masks | the vote map settings, and `plaque_vote_threshold`, `plaque_fp_min_votes`, `--min_size`, `--label_filter`, `label_filter_size` |
| vessels | `vessel_tile`, `vessel_stride`, `vessel_merge_mode`, `vessel_threshold`, the model URL (`models.vesselfm`) |

So a run that only changes the LANTERN vote threshold reuses the vote maps and
thresholds them again, without GPU inference. A distance transform is reused
only together with its mask. A prior dataset without `config/snakebids.yml`
has no segmentation or vessel outputs reused.

The brain mask that LANTERN and vesselFM use to skip empty tiles is not
compared: reuse the brain masks too (`brainmask`), or make them with the same
settings.
