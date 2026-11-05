# Registration Methods

Registration runs with a named method from `registration.presets` in
`config/snakebids.yml`. The name is also the `desc-` of the method's output
transforms, so results from different methods can sit side by side.

## Choosing a route

| Route | When | Needs |
|-------|------|-------|
| Register subjects to a packaged template such as ABAv3 | The default | CPU only |
| [Build a cohort template](cohort_template.md), then use it and its warps | Recommended when you can | A GPU, and at least 2 subjects |

A cohort template is not the default because of what it needs, but it is the
recommended route: subjects are registered to images of their own contrast, and the
hard step, from the cohort's shape to the atlas, is solved once on a high-SNR average.

## Options

| Option | Registration | Default |
|--------|--------------|---------|
| `--registration_method` | subject to `--template` | `greedy` |
| `--cohort-template-reg-method` | original template to cohort template | `antsMulti` |

The built-in methods are:

| Method | Backend | Description |
|--------|---------|-------------|
| `greedy` | greedy | NMI metric, the default |
| `greedyNoAffine` | greedy | no affine stage |
| `greedySharp` | greedy | less smoothing of the warp, to follow finer structure |
| `greedyNCC` | greedy | local cross-correlation metric, for images of similar contrast, such as a cohort template |
| `ants` | antsRegistration | SyN with the MI metric |
| `fireants` | FireANTs | MI loss (`fusedmi`); needs a GPU and the `gpu` environment |

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

Every method also has a multichannel variant, named `<method>Multi` (`greedyMulti`,
`antsMulti`, `fireantsMulti`, and so on), with an equally weighted metric term per
stain, so a structure that is only visible in one stain can still guide the
registration. Memory and run time grow with the number of stains.

For a cohort template registered to its original template, each cohort template
stain is compared against the same original template image. Compare these to their
base method with the labels in the report, since a stain with little in common with
the original template adds little.
`antsMulti` is the default here: it aligned best by eye when the methods above were
compared on one cohort template (`NkCreVATFX`, 5 subjects, 3 stains). Check the report
for your own data, as this is one comparison.

For subjects registered to a cohort template, each stain is compared against the
same stain of the template:

- With two or more stains shared with the template, `--registration_method` becomes
  its multichannel variant, using all of them. A template with two of a subject's
  three stains uses those two. `--registration_single_channel` turns this off.
- With one shared stain, registration uses that stain.
- With none, or a template without SPIM templates, registration uses the first
  stain in `--stains_for_reg` against the template's `anat` image.

Single-channel steps, such as the initial affine for brain masking, use the first
shared stain in `--stains_for_reg`, or else the first in the template. All stains
use that stain's brain mask. Warps from `--reuse_derivatives` are matched by method
name, so `desc-greedy` warps are not reused for a run that becomes `greedyMulti`.

## Comparing methods

These analysis levels run several methods and write one report with them side by
side. The report shows each method's atlas outlines, in all three planes, and can
flicker between the fixed and warped images. When the template has several stains,
you can switch between them, with every method shown on every stain:

```
spimquant /path/to/bids /path/to/output tune_reg \
    --registration_tune_methods greedy greedyMulti ants antsMulti \
    --participant-label 01 02

spimquant /path/to/bids /path/to/output cohort_template_tune_reg \
    --cohort-template-name lecanemab \
    --cohort-template-reg-tune-methods greedy greedyMulti ants antsMulti
```

The reports are `xfm/*_stain-<stain>_from-SPIM_to-<template>_desc-tune_regqc.html`
per subject, and
`tpl-<name>_stain-<stain>_from-<template>_to-<name>_desc-tune_regqc.html` for the
cohort template. Limit the subjects when tuning: every subject runs every method.
The metrics flag problems such as folding, but do not rank methods, since no metric
compares well across image contrasts. Pick by eye, then pass the method to
`--registration_method` or `--cohort-template-reg-method`.

The reports resample from the kept composite warps, so writing one again does not
repeat the registrations. The QC report of a single method also shows its affine
stage, which does need that registration's temporary outputs.

## Outputs

Each method writes `xfm/*_from-subject_to-<template>_desc-<method>_xfm.nii.gz` and
its inverse. `xfm/*_from-subject_to-<template>_xfm.nii.gz` and its inverse are links
to the warps selected for the subject, and are what the rest of the workflow uses.
