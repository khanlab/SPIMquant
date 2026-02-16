# Group Statistics

The group statistics workflow performs population-level analysis of segmentation metrics across multiple subjects.

## Purpose

Group-level analysis enables:
- Statistical comparison between experimental groups
- Identification of regions with significant differences
- Visualization of population-level pathology patterns
- Contrast-specific analysis and subgroup comparisons
- Creation of statistical parametric maps

The workflow integrates subject-level quantifications with participant metadata to perform region-wise statistical tests and generate comprehensive visualizations.

## Workflow Location

**File**: `spimquant/workflow/rules/groupstats.smk`

**Note**: This file already contains comprehensive docstrings for all rules describing their purpose and functionality.

## Workflow Stages

### 1. Subject Concatenation
Aggregates region-level statistics (segstats.tsv) or object-level data (parquet files) across all subjects, adding participant metadata.

### 2. Statistical Testing
Performs statistical comparisons between groups defined in participants.tsv. Tests are conducted independently for each brain region and metric.

### 3. Contrast-Specific Analysis
Creates subgroup-specific summaries and visualizations by filtering data based on experimental contrasts (e.g., control vs treatment).

### 4. Visualization
Generates heatmaps and volumetric statistical maps showing significant differences across brain regions.

### 5. Group-Level Density Maps
Creates population-level density maps by aggregating detected objects across subjects in template space.

## Key Rule Types

### Aggregation Rules
- Concatenate segmentation statistics across subjects
- Merge regionprops or colocalization data
- Join with participant metadata

### Statistical Testing Rules
- Group-wise comparisons per region
- Multiple comparison correction
- Effect size calculation

### Visualization Rules
- Heatmap generation for statistical results
- Volumetric statistical maps (NIfTI)
- Contrast-specific group averages

### Group Density Rules
- Aggregate point clouds across subjects
- Create population-level density maps
- Separate by contrast groups

## Statistical Methods

The workflow supports various statistical approaches:

### Group Comparison
- **T-tests**: Two-group comparisons
- **ANOVA**: Multi-group comparisons
- **Non-parametric tests**: For non-normal distributions

### Multiple Comparison Correction
- **False Discovery Rate (FDR)**: Controls expected proportion of false positives
- **Bonferroni**: Conservative family-wise error rate control
- **Benjamini-Hochberg**: Less conservative FDR procedure

### Effect Sizes
- **Cohen's d**: Standardized mean difference
- **Eta-squared**: Proportion of variance explained

## Configuration

### Contrast Definition

Define experimental groups in `participants.tsv`:

```tsv
participant_id    group    age    sex
sub-001          control   P60    F
sub-002          disease   P60    M
sub-003          control   P90    F
sub-004          disease   P90    M
```

Then in configuration:

```yaml
contrast_column: "group"
contrast_values: ["control", "disease"]
```

### Metrics for Analysis

```yaml
seg_metrics:
  - "count"
  - "density"
  - "fieldfrac"

coloc_seg_metrics:
  - "coloc_count"
  - "coloc_density"
```

## Data Types and Formats

### Inputs
- **Segmentation Statistics** (`.tsv`): Per-region metrics from all subjects
- **Region Properties** (`.parquet`): Object-level data from all subjects
- **Colocalization** (`.parquet`): Cross-channel associations
- **Participants Metadata** (`.tsv`): Group assignments and covariates

### Outputs
- **Group Statistics** (`.tsv`): Statistical test results per region
- **Heatmaps** (`.png`): Visual summaries of significant differences
- **Statistical Maps** (`.nii.gz`): Volumetric t-statistic or p-value maps
- **Group Average Maps** (`.nii.gz`): Mean metrics per contrast group
- **Concatenated Parquet** (`.parquet`): Combined object-level data

## Output Files

### Group Statistics TSV

Contains per-region statistics:

```tsv
region_id    region_name    metric         mean_control    mean_disease    t_stat    p_value    p_adj
123          Hippocampus    Abeta+density  42.3           85.7           -3.45      0.003      0.015
456          Cortex         Abeta+density  38.1           79.2           -2.98      0.008      0.032
```

### Statistical Maps

NIfTI files with brain regions colored by:
- T-statistics (for effect size)
- P-values (for significance)
- FDR-corrected p-values
- Effect sizes (Cohen's d)

## Common Workflows

### Two-Group Comparison

```bash
# Compare control vs disease groups
spimquant /path/to/data /path/to/output group \
  --contrast_column group \
  --contrast_values control disease
```

### Generate All Group Outputs

```yaml
# Request all group-level targets in config
targets_by_analysis_level:
  group:
    - "groupstats"        # Statistical comparisons
    - "groupavg"          # Group-specific averages  
    - "heatmaps"          # Visualization
    - "group_density"     # Population density maps
```

### Subgroup Analysis

```bash
# Analyze only young animals
spimquant /path/to/data /path/to/output group \
  --filter-participant age=P60
```

## Quality Control

### Visual Inspection

1. **Heatmaps**: Check for expected regional patterns
2. **Statistical maps**: Verify significance aligns with hypothesis
3. **Group density maps**: Compare pathology distribution between groups
4. **Effect sizes**: Assess biological meaningfulness

### Statistical Validation

- **Sample size**: Ensure adequate N per group (minimum 3-5)
- **Normality**: Check distribution of residuals
- **Outliers**: Identify and investigate extreme values
- **Multiple comparisons**: Verify appropriate correction applied

### Common Issues

**No significant differences**:
- Check sample size adequacy
- Verify group assignments in participants.tsv
- Assess metric variability
- Consider covariate adjustment

**Unexpected significant regions**:
- Check for batch effects
- Verify preprocessing consistency
- Inspect outliers
- Validate group definitions

**High false discovery rate**:
- Increase sample size
- Apply stricter FDR threshold
- Use more conservative correction method
- Validate preprocessing pipeline

## Interpretation Guidelines

### P-values
- **p < 0.05**: Nominally significant
- **p_adj < 0.05**: Significant after FDR correction
- **p_adj < 0.01**: Highly significant

### Effect Sizes
- **Cohen's d < 0.2**: Small effect
- **0.2 ≤ d < 0.8**: Medium effect
- **d ≥ 0.8**: Large effect

### Biological Relevance
Consider both statistical significance and effect size. Small p-values with tiny effect sizes may not be biologically meaningful.

## Performance Considerations

- **Memory**: Concatenating parquet files can require substantial RAM for large cohorts
- **Runtime**: Statistical tests are fast (<1 minute per metric)
- **Storage**: Group-level outputs are relatively small

## Dependencies

- **pandas**: Data manipulation and statistical operations
- **scipy**: Statistical testing functions
- **statsmodels**: Advanced statistical models
- **matplotlib/seaborn**: Visualization
- **nibabel**: NIfTI I/O for statistical maps

## Related Workflows

- **[Segmentation](segmentation.md)**: Produces per-subject statistics for group analysis
- **[Registration](templatereg.md)**: Enables spatial normalization for group comparisons
- **[Patches](patches.md)**: Can be used for qualitative assessment of group differences

## Example Analysis

### Complete Group Analysis Pipeline

1. **Process all subjects** to completion:
```bash
spimquant /path/to/data /path/to/output participant
```

2. **Define contrasts** in participants.tsv:
```tsv
participant_id    genotype    treatment
sub-001          WT          vehicle
sub-002          WT          drug
sub-003          KO          vehicle
sub-004          KO          drug
```

3. **Run group analysis**:
```bash
spimquant /path/to/data /path/to/output group \
  --contrast_column genotype \
  --contrast_values WT KO
```

4. **Examine outputs**:
- `groupstats.tsv`: Statistical results per region
- `groupstats.png`: Heatmap of significant regions
- `*_tstat.nii.gz`: Volumetric t-statistic map
- `*_groupavg.nii.gz`: Group-specific mean maps

### Multi-Factor Analysis

For complex designs (e.g., genotype × treatment), process contrasts separately:

```bash
# Genotype effect
spimquant /path/to/data /path/to/output group \
  --contrast_column genotype --contrast_values WT KO

# Treatment effect  
spimquant /path/to/data /path/to/output group \
  --contrast_column treatment --contrast_values vehicle drug
```

Then perform interaction analysis in external statistical software using the exported TSV files.
