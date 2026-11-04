# Registration Methods

Registration runs with a named method from `registration.presets` in
`config/snakebids.yml`. The name is also the `desc-` of the method's output
transforms, so results from different methods can sit side by side.

| Option | Registration | Default |
|--------|--------------|---------|
| `--registration_method` | subject to `--template` | `greedy` |
| `--cohort-template-reg-method` | original template to cohort template | `greedy` |

The built-in methods are:

| Method | Backend | Description |
|--------|---------|-------------|
| `greedy` | greedy | NMI metric, the default |
| `greedyNoAffine` | greedy | no affine stage |
| `greedySharp` | greedy | less smoothing of the warp, to follow finer structure |
| `greedyNCC` | greedy | local cross-correlation metric, for images of similar contrast |
| `greedySpim` | greedy | cohort template only: registers the original template's SPIM template with NCC |
| `ants` | antsRegistration | SyN with the MI metric |
| `antsSpim` | antsRegistration | cohort template only: registers the SPIM template with CC |
| `fireants` | FireANTs | MI loss; needs a GPU and the `gpu` environment |
| `fireantsSpim` | FireANTs | cohort template only: registers the SPIM template with CC |

To add a method, add a preset, starting from another one with `base:`. Names must
be alphanumeric:

```yaml
registration:
  presets:
    greedyFine:
      base: greedy
      deform_iters: 100x50x25
```

Greedy iterations are per level, with the last at full resolution, so `100x50` runs
at 2x downsampling and then at full resolution.

The `moving: spim` methods have the same contrast as the cohort template, but the
SPIM template may not be as well aligned to the atlases as the original template's
`anat`, so check the labels in the report.

## Comparing methods

These analysis levels run several methods and write one report with them side by
side. The report shows each method's atlas outlines, in all three planes, and can
flicker between the fixed and warped images:

```
spimquant /path/to/bids /path/to/output tune_reg \
    --registration_tune_methods greedy greedySharp ants fireants \
    --participant-label 01 02

spimquant /path/to/bids /path/to/output cohort_template_tune_reg \
    --cohort-template-name lecanemab \
    --cohort-template-reg-tune-methods greedy greedyNoAffine greedySpim
```

The reports are `xfm/*_stain-<stain>_from-SPIM_to-<template>_desc-tune_regqc.html`
per subject, and
`tpl-<name>_stain-<stain>_from-<template>_to-<name>_desc-tune_regqc.html` for the
cohort template. Limit the subjects when tuning: every subject runs every method.
The metrics flag problems such as folding, but do not rank methods, since no metric
compares well across image contrasts. Pick by eye, then pass the method to
`--registration_method` or `--cohort-template-reg-method`.

## Outputs

Each method writes `xfm/*_from-subject_to-<template>_desc-<method>_xfm.nii.gz` and
its inverse. `xfm/*_from-subject_to-<template>_xfm.nii.gz` and its inverse are links
to the warps selected for the subject, and are what the rest of the workflow uses.
