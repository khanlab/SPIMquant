# Landmark evaluation of registration (`--landmark_analysis`)

Measures registration accuracy with manual landmarks: target registration error
(TRE, µm) for **every registration method × every stain** SPIMquant can run, in one
participant-level run. It sits on top of the normal pipeline: the preprocessing
(downsampling, N4, brain mask) and the registrations are SPIMquant's own rules; the
landmark mode adds the method × stain loop, the point warping, the TRE and the report.

## What one run does

For each SPIM scan that has a landmarks file:

1. Builds the usual level-5 N4 + brain-masked image of **every stain** (one brain
   mask, from the default registration stain, shared by all stains). Scans may have
   different stains: each is registered with the stains it has, and the brain mask
   comes from a stain every scan has (the first such in `--stains_for_reg`, else the
   first shared stain).
2. Registers the scan to `--template` with **every method** in `--landmark_methods`
   (default: every single-channel `registration.presets` entry — greedy,
   greedyNoAffine, greedySharp, greedyNCC, ants, fireants, plus any preset added
   to the config) using **every stain** as the registration channel
   (`--landmark_stains`, default all stains; the `--stains_for_reg` priority is
   ignored).
3. With `--register_to_mri`: registers the MRI to SPIM once per stain (greedy,
   `reg_mri` config).
4. Computes TRE for every path:

| `path` | points | transform applied to the points | `inversion` |
|---|---|---|---|
| `spim_to_template` | SPIM → template | `from-<tpl>_to-subject` (shipped inverse) | `shipped` |
| `spim_to_template` | SPIM → template | fixed-point inverse of `from-subject_to-<tpl>` | `fixedpoint` |
| `template_to_spim` | template → SPIM | `from-subject_to-<tpl>` | `none` |
| `mri_to_spim` | MRI → SPIM | `from-SPIM_to-<mri>` | `none` |
| `spim_to_mri` | SPIM → MRI | `from-<mri>_to-SPIM` | `none` |
| `mri_to_template` | MRI → SPIM → template | both legs, same stain | `shipped` / `fixedpoint` |

Points move opposite to images: the transform that resamples B-space images into A
moves points from A to B. Landmarks are matched by name with any `Template_`,
`SPIM_`/`MRI_` prefix and `-<n>` suffix stripped; files are never modified.
Any number of landmarks can be used; each path uses the ones present on both ends.
Landmarks found on only one side (a misspelt name, or a point not placed in one set)
are left out of the TRE, so the run prints a warning listing them, and they are in
`group/landmark_name_check.tsv` and the report's *Landmark names* section.

## Usage

```bash
pixi run -e gpu spimquant <bids_dir> <output_dir> participant \
  --participant-label <labels> \
  --landmark_analysis \
  --landmarks_spim <manual_landmarks derivative>     # sub-*/[ses-*/]micr/*from-SPIM_landmarks.json
  --landmarks_template <folder of template points>   # markups in --template space
  [--register_to_mri --landmarks_mri <derivative>    # sub-*/[ses-*/]anat/*from-<suffix>_landmarks.json
   [--filter-mri ...] ]                               # session: from the landmark file
  [--landmark_methods greedy ants ...] [--landmark_stains YoPro Iba1 ...]
  [--landmark_keep_warps] [--landmark_regions regions.tsv]
```

- Use the `gpu` pixi environment and a GPU allocation if `fireants` is among the
  methods (default). Without a GPU, pass `--landmark_methods` without `fireants`.
- The MRI is chosen from the MRI landmark file: only the MRI scans matching the
  session (and any `acq-`/`run-`/`rec-`) in each subject's landmark file name are used,
  so other sessions (timepoints) are ignored without filtering them; the run prints
  the choice. If a subject has landmark files for several sessions, choose one with
  `--filter-mri session=<label>`. Several scans that still match (e.g. repeated runs
  of the session) are averaged on the grid of the first (sorted by name), so the MRI
  landmarks must be in that scan's space; put `run-<n>` in the landmark file name to
  use only that scan.
- Landmark files: 3D Slicer markups, `coordinateSystem` LPS (RAS is converted), mm.
  SPIM landmarks are on the level-5 grid of the scan (any stain).
- Subject→template composite warps are ~1.5 GB each at the 25 µm ABAv3 grid, so
  they are deleted after the TRE step unless `--landmark_keep_warps`.

## Outputs (`<output_dir>/landmarks/`)

- `sub-*/…/tre/*_tre.tsv` — one row per landmark × path × method × stain × inversion:
  `tre_um`, source/warped/target LPS coordinates, `inside_field` (False if a point
  fell outside a transform's grid, where ITK would use zero displacement), and the
  fixed-point residual.
- `group/landmarks_tre_allscans.tsv` — all scans.
- `group/landmarks_tre_summary.tsv` — n, median, IQR, mean, SD, max per
  path × method × stain × inversion (points outside a field excluded).
- `group/landmarks_tre_per_landmark.tsv` — the same per landmark.
- `sub-*/…/markups/` — where the landmarks land, as 3D Slicer markups (LPS mm),
  listed in `*_markups.tsv`:
  - `*_from-<A>_to-<B>_stain-<stain>_desc-<method>_landmarks.mrk.json`: the
    landmarks of space A carried to space B by that registration (orange); each
    point's description holds its TRE. SPIM→template and MRI→template use the
    fixed-point inverse.
  - `*_space-<B>_desc-target_landmarks.mrk.json`: the landmarks they are compared
    with (aqua).

  Open them in Slicer on the image of space B: the template (`tpl-<template>/…_anat`),
  the scan's `micr/*_level-5_desc-N4brain_SPIM.nii.gz` (SPIM), or the subject's
  `anat/*_desc-brain_<suffix>.nii.gz` (MRI).
- `group/report/report.md` (+ `figures/`, `*.tsv`, `markups/`) — landmark report,
  with a summary figure at the top (method, path, subject, region):
  0. failed registrations (median TRE > 2 mm on a scan, e.g. `greedyNoAffine`, which has
     no affine or centre initialisation and does not converge from the raw SPIM pose);
     listed and left out of the sections below;
  1. best method × stain per path (heatmaps);
  2. TRE per landmark (boxplots);
  3. TRE by region, with each path's best method and stain. Regions come from
     `--landmark_regions` (TSV with columns `region`, `pattern`: a landmark goes to the
     first region whose regular expression matches its name), or a built-in map for
     the cortex / hemisphere / ventricle / CC / striatum / thalamus / midbrain /
     cerebellum names (others are `other`);
  4. TRE by subject, for every path;
  5. landmark or registration? Per subject × landmark, the *consensus* (median
     warped position over all method × stain registrations) vs the landmark. A large
     consensus error with a small spread between registrations points at the landmark
     (placement or ambiguous anatomy); a large spread at the registration;
  6. consistency across subjects in template space (spread, bias of the mean from the
     template point, subject outliers: a subject's consensus point compared with the
     mean of the *other* subjects);
  7. left/right symmetry of paired landmarks (template and subjects, with a left vs
     right plot);
  8. MRI vs lightsheet routes to the template (table and per-landmark scatter);
  9. flagged points with screenshots (landmark on the SPIM image and on the template,
     with the consensus position). Flags only; landmarks are never edited.

  `report/markups/` holds the template points and every subject's consensus point in
  template space (labelled `<subject>_<landmark>`), to view all subjects at once.
- `xfm/` — the per-stain, per-method composite warps (inverse kept; forward only
  with `--landmark_keep_warps`) and MRI↔SPIM composites.

## Code

- `workflow/rules/landmarks.smk` — rules; the registration rules are the regular
  ones reused with `use rule … with:` and a `{regstain}` wildcard.
- `workflow/scripts/landmark_tre.py` — landmark I/O, point warping (SimpleITK,
  matches `antsApplyTransformsToPoints` to <0.05 µm), fixed-point inversion.
- `workflow/scripts/landmark_tre_subject.py`, `landmark_tre_summary.py`.
- `workflow/scripts/landmark_markups.py` — the per-scan Slicer markups.
- `workflow/scripts/landmark_report.py` — the report; also runs on its own
  (`python landmark_report.py --tsv … --template-landmarks … --out …`).
