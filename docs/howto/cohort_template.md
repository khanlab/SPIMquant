# Cohort templates

SPIMquant can build a study-specific (cohort) SPIM template from your subjects,
and prepare it for use as a `--template` in later SPIMquant runs. This is the
recommended way to register, but it is not the default, as it needs a GPU and at
least 2 subjects (see [Registration Methods](registration_methods.md#choosing-a-route)).

The workflow is:

1. Subjects are affine-registered to `--template` (e.g. ABAv3), and each stain is
   resampled to template space at `--cohort-template-level` (default:
   `--registration_level`). No deformable registration or earlier SPIMquant run is
   needed, and the template's shape comes only from the subjects.
2. A multivariate template is built from these images with
   [FireANTs](https://github.com/rohitrango/FireANTs) on GPU. All stains are built
   jointly, so every stain template is in the same space.
3. The cohort template has its own shape, so the original template is registered to the cohort template, with
   `--cohort-template-reg-method`, by default `antsMulti` (see
   [Registration Methods](registration_methods.md)). The original atlases, brain
   mask and anat are warped into cohort template space with this registration.
4. The `cohort_template_export` level writes a `template_description.json`, a
   SPIMquant `templates:` config entry, and optionally uploads everything to Zenodo.

## In one run: `--reg-cohort-template`

To register subjects to `--template` through a cohort template without the
separate `cohort_template` and `cohort_template_export` steps, add
`--reg-cohort-template` to a participant run:

```
spimquant /path/to/bids /path/to/output participant --reg-cohort-template
```

The cohort template is built in the same run, from all the subjects in it (steps 1
to 3 above, with the same `--cohort-template-*` options), and written to
`cohort_template/tpl-cohort/`. Each subject's warps to `--template` then chain its
affine to `--template`, its warp from the build, and the cohort template's
registration to `--template`, and are written as `desc-viacohort`. Unlike
`--template` with an exported cohort template, outputs stay in `--template` space
and use its atlases directly, so they can be compared with other runs registered
to it.

Since the template needs every subject, all subjects must be in the run (not one
`--participant_label` per job), and adding subjects later builds a new template.
To add subjects to a dataset without rebuilding, run with
`--reuse_derivatives <earlier output>`: subjects with `desc-viacohort` warps in it
reuse them, and the rest are registered with `--registration_method` (its
multichannel variant for a multi-stain template) to the cohort template in that
dataset, instead of a new one being built, so that every subject goes through
the same cohort template. See [Reusing Derivatives](reuse_derivatives.md).

The registration QC report of each subject (`xfm/*_desc-viacohort_regqc.html`)
shows its affine and its `desc-viacohort` warp, and the template's own
registration to `--template` has a report in `cohort_template/tpl-cohort/`.

`--template_crop` and cohort `--template`s are not supported with
`--reg-cohort-template`.

## Building the template

```
spimquant /path/to/bids /path/to/output cohort_template \
    --cohort-template-name lecanemab \
    --cohort-template-stains YoPro Iba1 Abeta \
    --cohort-template-n-gpus 2
```

This runs any registration that has not been run yet, so it can be run on a fresh
output directory or one with existing `participant` outputs.

Options:

- `--cohort-template-name`: the `tpl-` label of the new template (alphanumeric).
- `--cohort-template-stains`: stains to build from (default: all). The registration
  stain (e.g. YoPro) is always the first template channel, and is the image the
  original template is registered to. Stains with strong subject-specific signal
  (e.g. Abeta plaques) still contribute to the registration metric, so you may
  prefer to restrict this to anatomical stains.
- `--cohort-template-extra-inputs`: other SPIMquant output directories, from
  datasets registered to the same `--template`, to pool into the template. Subjects
  missing any of the stains are skipped.
- `--cohort-template-iterations`, `--cohort-template-n-gpus`: FireANTs settings.
  Other FireANTs parameters are in the `cohort_template:` section of
  `config/snakebids.yml`.

Outputs are in `<output>/cohort_template/`:

```
tpl-lecanemab/
├── tpl-lecanemab_stain-YoPro_SPIM.nii.gz          # one per stain, normalized to [0,1]
├── tpl-lecanemab_desc-rawintensity_stain-YoPro_SPIM.nii.gz   # with --cohort-template-warp-raw
├── tpl-lecanemab_anat.nii.gz                      # original anat, warped
├── tpl-lecanemab_desc-brain_mask.nii.gz           # original mask, warped
├── tpl-lecanemab_seg-<seg>_dseg.nii.gz/.tsv       # original atlases, warped
├── tpl-lecanemab_from-ABAv3_to-lecanemab_desc-greedy_xfm.nii.gz
├── tpl-lecanemab_from-lecanemab_to-ABAv3_desc-greedy_xfm.nii.gz
└── tpl-lecanemab_stain-YoPro_from-ABAv3_to-lecanemab_desc-greedy_regqc.html
tpl-ABAv3/
└── tpl-ABAv3_desc-lecanemab_stain-YoPro_SPIM.nii.gz   # cohort template in ABAv3 space
```

Check the `regqc.html` report to confirm the atlases line up with the cohort template.

Template building clips each subject image to robust percentiles and rescales it
to [0,1], and FireANTs sharpens and min-max normalizes the template, so the stain
templates do not keep the original intensities. The `desc-rawintensity` templates,
made with `--cohort-template-warp-raw`, do: each is the voxelwise mean of the subjects' un-normalized images (the
`space-<template>_desc-affine` SPIM images, which have no bias field correction or other
intensity changes), warped with the warps from the final template iteration. They
are in the same space as the stain templates but are not sharpened, so they look
slightly smoother. Use them when absolute intensity matters, e.g. with uniform
laser power settings across subjects; the stain templates remain the ones used for
registration.

FireANTs is only installed in the `gpu` pixi environment, so run this level from
it, e.g. `pixi run -e gpu spimquant ...`. It uses the fused CUDA ops if they are
installed, falling back to the plain implementation otherwise.

Shape averaging (`cohort_template.fireants.shape_avg`) is off by default: it adds
a full-resolution registration per iteration, runs out of memory at the 25um
ABAv3 resolution on a 40GB GPU, and in testing gave a worse subject fit than the
plain average. To lower memory use and run time, resample the subject images
before template building with e.g. `--cohort-template-resample-mm 0.05`, useful
when using a finer `--cohort-template-level`. On an 18-subject cohort, a 50um
template built in ~7 min vs ~50 min at 25um.

The fused CUDA ops (`fireants-fused-ops`, from the `akhanf/khanlab` channel) are
installed in the gpu environment and reduce FireANTs' GPU memory use.

## Preparing the template for SPIMquant

```
spimquant /path/to/bids /path/to/output cohort_template_export \
    --cohort-template-name lecanemab \
    --cohort-template-stains YoPro Iba1 Abeta
```

Use the same `cohort_template` options as when building. This writes
`tpl-lecanemab/template_description.json`, `tpl-lecanemab_buildinfo.tar.gz` (image
list, FireANTs config and logs), and `tpl-lecanemab_spimquant.yml`, a SPIMquant
config file that adds the template and selects it.

To upload to Zenodo, create a draft upload on Zenodo and pass its deposition ID,
with a Zenodo access token in the environment:

```
export SNAKEMAKE_STORAGE_ZENODO_ACCESS_TOKEN=<token>
spimquant /path/to/bids /path/to/output cohort_template_export \
    --cohort-template-name lecanemab \
    --cohort-template-stains YoPro Iba1 Abeta \
    --cohort-template-zenodo-deposition 12345678
```

The config file then points at `https://zenodo.org/records/<id>/files/...` URLs,
which resolve once you publish the deposition. Without a deposition, it points at
the absolute paths of the exported files.

## Using the template

Pass the config file to any SPIMquant run, in place of `--template`:

```
spimquant /path/to/bids /path/to/output participant \
    --configfile tpl-lecanemab_spimquant.yml
```

The registration stain must be one of the cohort template's stains (see
`--stains_for_reg`).

Packaged templates such as ABAv3 are registered to through their `anat` image. Only
cohort templates define stain templates. Subjects are registered directly to the
ones they have a stain for, using all of them when there are two or more (see
[Registration Methods](registration_methods.md)).

Subjects that were used to build the template are not registered again. Their
warps to the template are composed from the warps saved by the build, and from
their `desc-affine` warps to the original template in the dataset the template was
built in, and are written as `desc-cohort`. Other subjects are registered with
`--registration_method`, or its multichannel variant. This needs both of those to still be in place at
the paths in `tpl-lecanemab_spimquant.yml`, and a template built with a version of
SPIMquant that saves the inverse warps (`*_invwarp.nii.gz` in the build directory);
otherwise the subject is registered as usual.

To make the template a built-in `--template` choice, open a PR that moves its
`templates:` entry into `spimquant/config/snakebids.yml` (with the Zenodo URLs, or
paths relative to `spimquant/resources/` after copying `tpl-lecanemab` there), and
adds the name to the `--template` choices.
