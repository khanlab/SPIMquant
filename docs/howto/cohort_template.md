# Cohort templates

SPIMquant can build a study-specific (cohort) SPIM template from your subjects,
and prepare it for use as a `--template` in later SPIMquant runs.

The workflow is:

1. Subjects are registered to `--template` (e.g. ABAv3) as in the `participant`
   level, and each stain is resampled to template space at `--cohort-template-level`
   (default: `--registration_level`).
2. A multivariate template is built from these images with
   [FireANTs](https://github.com/rohitrango/FireANTs) on GPU. All stains are built
   jointly, so every stain template is in the same space.
3. The cohort template drifts away from the original template's space during
   building, so the original template's anat is registered (affine + deformable,
   greedy) to the cohort template. The original atlases, brain mask and anat are
   warped into cohort template space with this registration.
4. The `cohort_template_export` level writes a `template_description.json`, a
   SPIMquant `templates:` config entry, and optionally uploads everything to Zenodo.

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
├── tpl-lecanemab_stain-YoPro_SPIM.nii.gz          # one per stain
├── tpl-lecanemab_anat.nii.gz                      # original anat, warped
├── tpl-lecanemab_desc-brain_mask.nii.gz           # original mask, warped
├── tpl-lecanemab_seg-<seg>_dseg.nii.gz/.tsv       # original atlases, warped
├── tpl-lecanemab_from-ABAv3_to-lecanemab_xfm.nii.gz
├── tpl-lecanemab_from-lecanemab_to-ABAv3_xfm.nii.gz
└── tpl-lecanemab_stain-YoPro_from-ABAv3_to-lecanemab_regqc.html
tpl-ABAv3/
└── tpl-ABAv3_desc-lecanemab_stain-YoPro_SPIM.nii.gz   # cohort template in ABAv3 space
```

Check the `regqc.html` report to confirm the atlases line up with the cohort template.

The FireANTs step runs in its own conda environment (`workflow/envs/fireants.yaml`),
and uses the fused CUDA ops if they are installed, falling back to the plain
implementation otherwise.

## Preparing the template for SPIMquant

```
spimquant /path/to/bids /path/to/output cohort_template_export \
    --cohort-template-name lecanemab \
    --cohort-template-stains YoPro Iba1 Abeta
```

Use the same `cohort_template` options as when building. This writes
`tpl-lecanemab/template_description.json`, `tpl-lecanemab_buildinfo.tar.gz` (image
list, FireANTs config and logs), and `tpl-lecanemab_spimquant.yml`, a `templates:`
entry to add to `spimquant/config/snakebids.yml` (also add the name to the
`--template` choices).

To upload to Zenodo, create a draft upload on Zenodo and pass its deposition ID,
with a Zenodo access token in the environment:

```
export SNAKEMAKE_STORAGE_ZENODO_ACCESS_TOKEN=<token>
spimquant /path/to/bids /path/to/output cohort_template_export \
    --cohort-template-name lecanemab \
    --cohort-template-stains YoPro Iba1 Abeta \
    --cohort-template-zenodo-deposition 12345678
```

The config entry then points at `https://zenodo.org/records/<id>/files/...` URLs,
which resolve once you publish the deposition. Without a deposition, it points at
`tpl-lecanemab/...` paths relative to `spimquant/resources/`, where you can copy the
`tpl-lecanemab` directory instead.
