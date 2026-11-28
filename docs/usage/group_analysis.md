# Group analysis

The `group` level combines the participant outputs in an output directory
and, optionally, compares groups region by region. Run it after
`participant`, on the same output directory:

```bash
spimquant bids/ out/ group --cores all
```

Pass the same options as to the participant run (levels, atlases, methods,
`--no_vessels`...): the group level checks that the participant outputs are
up to date, and with different options it plans, and runs, a different
participant level.

Without a model this only writes the combined tables. All group outputs go
to `out/group/<label>/` (`--group-stats-label`, default `1`), so analyses with
different models or subjects do not overwrite each other.

## participants.tsv

The columns of `participants.tsv` in the BIDS directory are joined onto every
row of the group tables and are the variables of the model:

```tsv
participant_id	treatment	genotype	sex	age
sub-01	vehicle	WT	M	12
sub-02	drug	WT	F	13
sub-03	vehicle	KO	M	11
```

## Combined tables

| File | Rows |
|---|---|
| `seg-{atlas}_from-{template}_allsubjects.tsv` | subject × region: every column of `mergedsegstats.tsv` plus the participant columns |
| `stain-{stain}_..._allsubjects_regionpropstats.tsv` | every object of one stain in every subject, with its region and participant columns |
| `space-{template}_regionprops.parquet` | every object of every stain in every subject, with distances |
| `space-{template}_coloc.parquet` | every colocalized pair in every subject |

These are the tables to load for analyses of your own (R, pandas, Prism).
Each has a JSON sidecar describing its columns.

## Group QC report

Every group run also writes `group/qc/groupqc.html`: the registration and
segmentation QC of all subjects on one page, with a montage, an animated
flip through the subjects and outlier flags (`group/qc/groupqc.tsv`). Check
it before the statistics; see [Group QC report](../howto/group_qc.md).

## Region-wise statistics

```bash
spimquant bids/ out/ group \
  --group-stats-model "metric ~ C(treatment) + age + C(sex)" \
  --group-stats-pairwise treatment
```

| Option | |
|---|---|
| `--group-stats-model` | A [patsy](https://patsy.readthedocs.io/) formula with `metric` as the response. It is fitted by ordinary least squares, separately for every region and every metric (`Abeta+fieldfrac`, `Iba1+density`, ...). |
| `--group-stats-pairwise` | Factor(s) whose levels are compared pairwise. Repeat for several factors. |
| `--group-stats-within` | Factor(s) to stratify by: each contrast is computed within every combination of their levels (e.g. `genotype sex` gives WT/M, WT/F, KO/M, KO/F). |
| `--group-stats-where` | A [pandas query](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) that selects the subjects used for fitting and contrasts, e.g. `"treatment in ['vehicle', 'drug'] and qc_pass == 1"`. |
| `--group-stats-label` | Name of the output folder under `group/`. |

For each contrast, e.g. `contrast-treatment+vehiclevsdrug` (or
`contrast-treatment+vehiclevsdrug+genotype-WT+sex-M` with strata), it writes:

| File | Content |
|---|---|
| `*_contrast-{contrast}_groupstats.tsv` | One row per region: group sizes, and for each metric the two model-predicted group means, t-statistic, p-value and Cohen's d. |
| `*_contrast-{contrast}_groupstats.png` | Heatmap of the t-statistics, regions × metrics. |
| `*_contrast-{contrast}_metric-{metric}_{tstat,pval,cohensd}.nii.gz` | The same values painted into the template atlas, one image per metric and statistic. |

How a contrast is computed: the model is fitted on the subjects passing
`--group-stats-where`. The means of the two levels are predicted with
continuous covariates at their mean, other categorical factors at their most
frequent level, and strata at the requested level. The difference of the two
predictions is tested against the model's covariance (t-statistic, p-value).
Cohen's d uses the pooled standard deviation of the two groups.

!!! warning "Multiple comparisons"
    p-values are not corrected. With hundreds of regions and several metrics,
    correct them yourself, e.g. with `statsmodels.stats.multitest.multipletests`
    in Python or `p.adjust()` in R, on the `groupstats.tsv` columns.

## Examples

Genotype within each sex:

```bash
spimquant bids/ out/ group --group-stats-label genotype_by_sex \
  --group-stats-model "metric ~ C(genotype) * C(sex) + age" \
  --group-stats-pairwise genotype --group-stats-within sex
```

Two treatment arms, excluding subjects that failed QC:

```bash
spimquant bids/ out/ group --group-stats-label lecanemab \
  --group-stats-model "metric ~ C(treatment) * C(sex) + age" \
  --group-stats-pairwise treatment \
  --group-stats-where "treatment in ['PBS', 'Lecanemab'] and qc_pass == 1"
```
