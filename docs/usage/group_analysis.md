# Group-Level Statistical Analysis

SPIMquant supports group-level statistical analysis to compare pathology metrics across experimental groups using a formula-based OLS model (via [statsmodels](https://www.statsmodels.org/) and [patsy](https://patsy.readthedocs.io/)).

## Overview

Group analysis performs statistical comparisons of:

- **Field fraction**: Proportion of pathology signal per atlas region
- **Density**: Cell/particle count per unit volume
- **Count / Volume / Voxel counts**: Other quantification metrics
- **Colocalization metrics**: Overlap ratio, distance, density, count

Results are provided as:

- Merged per-subject ROI tables for export to external stats tools
- Statistical tables per pairwise contrast (t-statistics, p-values, Cohen's d, group means)
- Heatmap visualizations
- 3D volumetric stat maps

## Prerequisites

Before running group analysis:

1. **Complete participant-level processing** for all subjects
2. **Create participants.tsv** with group assignments

## Creating participants.tsv

Create a tab-separated file in your BIDS root directory with participant metadata:

```tsv
participant_id	treatment	genotype	sex	age
sub-01	vehicle	WT	M	12
sub-02	vehicle	WT	F	13
sub-03	drug	WT	M	11
sub-04	drug	WT	F	12
sub-05	vehicle	KO	M	12
sub-06	vehicle	KO	F	14
sub-07	drug	KO	M	11
sub-08	drug	KO	F	12
```

Required columns:

- `participant_id`: Subject identifiers matching BIDS subject IDs (e.g. `sub-01`)
- At least one grouping column (e.g., `treatment`, `genotype`)

## Running Group Analysis

### Minimal Example

Fit a simple two-group model and compute all pairwise contrasts for `treatment`:

```bash
pixi run spimquant /bids /output group \
  --group-stats-model "metric ~ C(treatment)" \
  --group-stats-pairwise treatment \
  --cores all
```

### Including Covariates

Add continuous or categorical covariates to the model:

```bash
pixi run spimquant /bids /output group \
  --group-stats-model "metric ~ C(treatment) + age + C(sex)" \
  --group-stats-pairwise treatment \
  --cores all
```

### Interaction Effects

Fit a full factorial model with interaction terms:

```bash
pixi run spimquant /bids /output group \
  --group-stats-model "metric ~ C(treatment) * C(genotype) * C(sex) + age" \
  --group-stats-pairwise treatment \
  --group-stats-pairwise genotype \
  --cores all
```

Multiple `--group-stats-pairwise` flags enumerate pairwise contrasts for each factor independently.

### Stratified Contrasts

Use `--group-stats-within` to compute treatment contrasts *within* each combination of stratifying factors (e.g., separately for each genotype × sex combination):

```bash
pixi run spimquant /bids /output group \
  --group-stats-model "metric ~ C(treatment) * C(genotype) * C(sex) + age" \
  --group-stats-pairwise treatment \
  --group-stats-within genotype sex \
  --cores all
```

This produces one set of output files per stratum combination (e.g., `WT`/`M`, `WT`/`F`, `KO`/`M`, `KO`/`F`).

### Filtering Subjects

Use `--group-stats-where` to restrict which subjects are included in model fitting and contrasts.
The expression is a [pandas query string](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html):

```bash
# Exclude subjects not meeting QC criteria
pixi run spimquant /bids /output group \
  --group-stats-model "metric ~ C(treatment) + age" \
  --group-stats-pairwise treatment \
  --group-stats-where "treatment in ['vehicle', 'drug'] and qc_pass == 1" \
  --cores all
```

Subjects that do not match the expression are excluded from both model fitting and contrast enumeration.

### Labelling Analysis Runs

Use `--group-stats-label` to assign a name to each analysis run. Outputs are placed under `<output_dir>/group/<label>/`, so different analyses (different models, different filters, etc.) never overwrite each other:

```bash
# First analysis
pixi run spimquant /bids /output group \
  --group-stats-label treatment_only \
  --group-stats-model "metric ~ C(treatment)" \
  --group-stats-pairwise treatment \
  --cores all

# Second analysis with covariates
pixi run spimquant /bids /output group \
  --group-stats-label treatment_with_covariates \
  --group-stats-model "metric ~ C(treatment) + age + C(sex)" \
  --group-stats-pairwise treatment \
  --cores all
```

Default label is `1` if not specified.

## Group Analysis Options Summary

| Option | Description |
|--------|-------------|
| `--group-stats-model` | Patsy/statsmodels formula with `metric` as response placeholder, e.g. `"metric ~ C(treatment) + age"` |
| `--group-stats-pairwise` | Factor(s) for which all pairwise comparisons are computed. Can be repeated for multiple factors. |
| `--group-stats-within` | Stratifying factor(s) — contrasts are computed separately within each level combination. |
| `--group-stats-where` | Pandas query expression to filter subjects before model fitting. |
| `--group-stats-label` | Label for this run; outputs go under `group/<label>/` (default: `1`). |

## Output Files

All group outputs are placed under `<output_dir>/group/<label>/`.

### Merged Subject Table (always produced)

**`*_allsubjects.tsv`**: All subject-level ROI statistics concatenated into a single file, joined with participant metadata from `participants.tsv`. Produced whenever `analysis_level=group`, regardless of whether any contrasts are specified.

```tsv
index	name	Abeta+fieldfrac	Abeta+density	...	participant_id	treatment	sex	age
0	Isocortex	0.023	1452	...	sub-01	vehicle	M	12
0	Isocortex	0.031	1893	...	sub-03	drug	M	11
...
```

Useful for exporting to R, Python (pandas/seaborn), or GraphPad Prism for custom analysis.

An accompanying **`*_allsubjects.json`** sidecar describes each column.

### Statistical Results (per pairwise contrast)

**`*_contrast-<label>_groupstats.tsv`**: Per-region statistics for one pairwise contrast.

The contrast label encodes the factor, levels, and any strata — for example:
- `contrast-treatment+vehiclevsdrug` — vehicle vs. drug, no strata
- `contrast-treatment+vehiclevsdrug+genotype-WT+sex-M` — same contrast within WT males

```tsv
index	name	n_vehicle	n_drug	Abeta+fieldfrac_mean_vehicle	Abeta+fieldfrac_mean_drug	Abeta+fieldfrac_tstat	Abeta+fieldfrac_pval	Abeta+fieldfrac_cohensd	...
0	Isocortex	10	10	0.023	0.031	-2.14	0.042	-0.95	...
1	Hippocampus	10	10	0.015	0.024	-1.88	0.071	-0.84	...
```

Columns per metric:
- `<metric>_mean_<levelA>` / `<metric>_mean_<levelB>`: Group means
- `<metric>_tstat`: OLS-derived t-statistic for the marginal mean contrast
- `<metric>_pval`: Corresponding p-value (uncorrected)
- `<metric>_cohensd`: Cohen's d effect size

An accompanying **`*_groupstats.json`** sidecar describes each column.

### Visualizations

**`*_contrast-<label>_groupstats.png`**: Heatmap of t-statistics across all brain regions and metrics for one contrast.

### 3D Maps

**`*_contrast-<label>_groupstats.nii`**: Volumetric map of t-statistics registered to the template space.

View in FSLeyes, ITK-SNAP, or 3D Slicer.

## Statistical Methods

SPIMquant fits a single **OLS (Ordinary Least Squares)** model per brain region and metric using the formula you supply:

1. The user formula (e.g. `metric ~ C(treatment) + age`) is evaluated for each region × metric combination using patsy for design matrix construction.
2. The model is fit on all rows passing the `--group-stats-where` filter.
3. **Marginal means** for the two contrast levels are predicted at reference covariate values (continuous covariates at their mean, categorical at their mode; strata factors fixed at the requested level).
4. The **contrast vector** (levelA − levelB) is evaluated against the model's covariance matrix to produce the t-statistic and p-value.
5. **Cohen's d** is computed from pooled standard deviation.

!!! note "Multiple comparison correction"
    P-values in the output are uncorrected. Apply FDR or Bonferroni correction in your downstream analysis (e.g. using `statsmodels.stats.multitest.multipletests` in Python or `p.adjust()` in R).

## Example Workflows

### Drug Treatment Study

```bash
# 1. Process all subjects (participant level)
pixi run spimquant /bids /output participant --cores all

# 2. Run group analysis comparing drug vs. vehicle
pixi run spimquant /bids /output group \
  --group-stats-label drug_study \
  --group-stats-model "metric ~ C(treatment) + age + C(sex)" \
  --group-stats-pairwise treatment \
  --cores all
```

### Genotype Comparison with Sex Stratification

```bash
pixi run spimquant /bids /output group \
  --group-stats-label genotype_by_sex \
  --group-stats-model "metric ~ C(genotype) * C(sex) + age" \
  --group-stats-pairwise genotype \
  --group-stats-within sex \
  --cores all
```

This produces separate statistical maps for male and female cohorts.

### Multi-Factor Study with Subject Exclusion

```bash
pixi run spimquant /bids /output group \
  --group-stats-label lecanemab_qc \
  --group-stats-model "metric ~ C(treatment) * C(sex) + age" \
  --group-stats-pairwise treatment \
  --group-stats-where "treatment in ['PBS', 'Lecanemab'] and qc_pass == 1" \
  --cores all
```

### Merged Table Only (No Statistics)

To obtain only the merged subject table without running any statistical contrasts, simply omit `--group-stats-pairwise`:

```bash
pixi run spimquant /bids /output group --cores all
```

This always runs when `analysis_level=group` and produces `*_allsubjects.tsv` for all subjects.

## Interpreting Results

### Statistical Significance

- Interpret p-values in the context of the number of regions and metrics tested
- Apply appropriate multiple comparison correction for your use case
- Consider effect sizes (Cohen's d) alongside p-values

### Effect Sizes (Cohen's d)

| d | Interpretation |
|---|----------------|
| < 0.2 | Small |
| 0.2 – 0.8 | Medium |
| ≥ 0.8 | Large |

### Heatmap Color Scale

- **Positive t-statistic** (warm colors): levelA > levelB
- **Negative t-statistic** (cool colors): levelA < levelB
- Intensity reflects the magnitude of the difference

## Troubleshooting

### Missing participants.tsv

Error: `participants.tsv not found`

- Create `participants.tsv` in the BIDS root directory
- Ensure tab-separated format

### NaN in tstat / pval columns

Check the log file at `logs/group_stats/<entities>_groupstats.log`. Common causes:
- Fewer than 2 subjects per group in a region (logged as a skip message)
- A column referenced in the formula is missing from the data
- A contrast level is absent after `--group-stats-where` filtering

### Contrast Level Not Found

Error: `Contrast level 'X' ... is not present in the filtered data`

- Verify the spelling of level names in `participants.tsv`
- Check that `--group-stats-where` does not exclude subjects with that level

### Group Size Warnings

Minimum 2 subjects per group are required to fit the model. Larger groups increase statistical power.

## Next Steps

- [CLI Reference](cli.md): Full group analysis option reference
- [Visualization](../howto/visualization.md): Advanced visualization techniques
- [FAQ](../faq.md): Common questions about group analysis
