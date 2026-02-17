# Group-Level Statistical Analysis

SPIMquant supports group-level statistical analysis to compare pathology metrics across experimental groups.

## Overview

Group analysis performs statistical comparisons of:

- **Field fraction**: Proportion of pathology signal
- **Density**: Count-based metrics
- **Volume**: Volumetric measurements
- **Custom metrics**: User-defined quantifications

Results are provided as:

- Statistical tables (t-statistics, p-values, effect sizes)
- Heatmap visualizations
- 3D volumetric maps
- Group-averaged statistics

## Prerequisites

Before running group analysis:

1. **Complete participant-level processing** for all subjects
2. **Create participants.tsv** with group assignments
3. Define contrast column and values

## Creating participants.tsv

Create a TSV file in your BIDS directory with participant metadata:

```tsv
participant_id	treatment	age	sex	genotype
sub-01	control	12	M	WT
sub-02	control	13	F	WT
sub-03	drug	11	M	WT
sub-04	drug	12	F	WT
sub-05	control	12	M	KO
sub-06	drug	11	M	KO
```

Required columns:

- `participant_id`: Subject identifiers (must match BIDS subject IDs)
- At least one grouping column (e.g., `treatment`, `genotype`)

## Running Group Analysis

### Basic Group Comparison

Compare two groups:

```bash
pixi run spimquant /bids /output group \
  --contrast_column treatment \
  --contrast_values control drug \
  --cores all
```

This compares:

- **Group 1**: `treatment == control`
- **Group 2**: `treatment == drug`

### Multiple Contrasts

<!-- TODO: Document multiple contrast support -->

```bash
# Compare multiple grouping variables
pixi run spimquant /bids /output group \
  --contrast_column genotype \
  --contrast_values WT KO \
  --cores all
```

## Output Files

Group analysis generates several output files:

### Statistical Results

**`*_groupstats.tsv`**: Statistical test results by region

```tsv
region	t_statistic	p_value	effect_size	mean_control	mean_drug
ctx-lh-frontal	3.45	0.023	0.82	0.145	0.234
ctx-lh-temporal	2.11	0.058	0.51	0.098	0.143
...
```

Columns:

- `region`: Brain region name
- `t_statistic`: t-test statistic
- `p_value`: Uncorrected p-value
- `effect_size`: Cohen's d effect size
- `mean_<group>`: Group mean values

### Visualizations

**`*_groupstats.png`**: Heatmap of statistical results

Shows color-coded significance across brain regions.

### 3D Maps

**`*_groupstats.nii`**: Volumetric maps of statistics

Contains:

- t-statistic maps
- p-value maps
- Effect size maps

View in neuroimaging software (FSLeyes, ITK-SNAP, 3D Slicer).

### Group Averages

**`*_groupavgsegstats.tsv`**: Average statistics per group

```tsv
region	group	mean_fieldfrac	std_fieldfrac	n_subjects
ctx-lh-frontal	control	0.145	0.023	10
ctx-lh-frontal	drug	0.234	0.034	10
...
```

**`*_groupavg.nii.gz`**: Group-averaged volumetric maps

Separate volumes for each group showing average pathology distribution.

## Statistical Methods

<!-- TODO: Add details on statistical methods -->

SPIMquant uses:

- **t-tests**: For two-group comparisons
- **Multiple comparison correction**: Optional FDR/Bonferroni
- **Effect size**: Cohen's d for interpretability

## Filtering and Quality Control

### Excluding Subjects

<!-- TODO: Document subject exclusion -->

Exclude specific subjects from analysis:

```bash
pixi run spimquant ... group \
  --exclude_subjects 05 07 \
  --contrast_column treatment \
  --contrast_values control drug
```

### Region Filtering

<!-- TODO: Document region filtering -->

Focus analysis on specific brain regions.

## Visualization Options

<!-- TODO: Add visualization customization -->

### Customize Heatmaps

Control heatmap appearance:

```bash
# TODO: Add heatmap customization options
```

### Export for External Tools

Results can be imported into:

- **R**: For custom statistical analysis
- **Python**: Using pandas/seaborn
- **GraphPad Prism**: For publication figures

## Advanced Analysis

### Covariates

<!-- TODO: Document covariate support -->

Include covariates in statistical models:

```bash
# TODO: Add covariate examples
```

### Multiple Testing Correction

<!-- TODO: Document multiple testing correction -->

Apply FDR or Bonferroni correction:

```bash
# TODO: Add correction examples
```

### Custom Contrasts

<!-- TODO: Document custom contrasts -->

Define complex contrasts:

```bash
# TODO: Add custom contrast examples
```

## Example Workflows

### Drug Treatment Study

```bash
# 1. Process all subjects
pixi run spimquant /bids /output participant --cores all

# 2. Run group analysis
pixi run spimquant /bids /output group \
  --contrast_column treatment \
  --contrast_values vehicle drug \
  --cores all
```

### Genotype Comparison

```bash
# Compare wildtype vs knockout
pixi run spimquant /bids /output group \
  --contrast_column genotype \
  --contrast_values WT KO \
  --cores all
```

### Multi-Factor Design

<!-- TODO: Add multi-factor examples -->

```bash
# TODO: Add interaction effects
```

## Interpreting Results

### Statistical Significance

- **p < 0.05**: Traditionally significant
- **p < 0.01**: Highly significant
- Consider effect sizes alongside p-values

### Effect Sizes

Cohen's d interpretation:

- **d < 0.2**: Small effect
- **0.2 ≤ d < 0.8**: Medium effect
- **d ≥ 0.8**: Large effect

### Visualization

<!-- TODO: Add interpretation guide -->

- Red/warm colors: Higher in group 2
- Blue/cool colors: Higher in group 1
- Intensity: Magnitude of difference

## Troubleshooting

### Missing participants.tsv

Error: `participants.tsv not found`

- Create participants.tsv in BIDS root directory
- Ensure TSV format (tab-separated)

### Group Size Warnings

Warnings about small group sizes:

- Minimum 3 subjects per group recommended
- Larger groups increase statistical power

### Incomplete Data

Some subjects missing outputs:

- Ensure all subjects completed participant-level
- Check for failed jobs in participant-level logs

## Next Steps

- [Visualization](../howto/visualization.md): Advanced visualization techniques
- [Examples](../examples/workflows.md): Complete analysis examples
- [FAQ](../faq.md): Common questions about group analysis