# Documentation Scripts

This directory contains scripts for generating and maintaining SPIMquant documentation.

## Scripts

### `generate_dag_diagrams.py`

Generates modular Mermaid diagrams from the SPIMquant Snakemake workflow.

**Purpose:**
- Generate visual representations of the workflow DAG
- Break down the complex workflow into digestible, stage-specific diagrams
- Enable automated regeneration of diagrams as the workflow evolves

**Usage:**

```bash
# Basic usage (uses tests/bids_ds by default)
python3 generate_dag_diagrams.py

# With custom BIDS dataset
python3 generate_dag_diagrams.py --bids-dir /path/to/bids/dataset

# With different template
python3 generate_dag_diagrams.py --template gubra

# Full options
python3 generate_dag_diagrams.py --bids-dir /path/to/bids --output-dir /tmp/out --template ABAv3
```

**Options:**
- `--bids-dir PATH` - Path to BIDS dataset (default: `tests/bids_ds`)
- `--output-dir PATH` - Path to workflow output directory (default: `/tmp/spimquant_output`)
- `--template STR` - Template to use for registration (default: `ABAv3`)

**Output:**
- `docs/figures/rulegraph_full.mermaid` - Complete workflow rulegraph
- `docs/figures/dag_*.mermaid` - Individual stage diagrams (11 stages)

**Workflow Stages:**

The script categorizes workflow rules into these functional stages:

1. **01_import** - Import and setup templates, masks, and atlases
2. **02_preprocessing** - Image preprocessing and downsampling
3. **03_masking** - Brain masking using atropos segmentation
4. **04_correction** - Intensity correction and normalization
5. **05_registration** - Template registration (affine and deformable)
6. **06_transform** - Apply transformations to images and atlases
7. **07_segmentation** - Segmentation of pathology (threshold, multi-otsu)
8. **08_quantification** - Region properties and quantification
9. **09_statistics** - Statistical analysis and atlas mapping
10. **10_qc** - Quality control and reporting
11. **11_patches** - Extract image patches for analysis

**How It Works:**

1. **Generate Rulegraph**: Runs `spimquant.run` with `--rulegraph mermaid-js` flag to get the complete workflow graph
2. **Parse Mermaid**: Extracts nodes (rules) and edges (dependencies) from the Mermaid flowchart
3. **Classify Nodes**: Uses regex patterns to classify each rule into a functional stage
4. **Create Subgraphs**: For each stage, generates a diagram showing:
   - Primary rules for that stage (highlighted in green)
   - Dependent rules from other stages (shown in light blue)
   - All relevant edges connecting the rules
5. **Save Diagrams**: Writes Mermaid files to `docs/figures/`

**Dependencies:**

Required Python packages:
- `snakemake` (>=9.0)
- `snakebids` (>=0.14.0)

**Adding New Stages:**

To add or modify workflow stages, edit the `WORKFLOW_STAGES` dictionary in the script:

```python
WORKFLOW_STAGES = {
    "stage_name": {
        "description": "Stage description",
        "patterns": [
            r"^rule_name_pattern$",
            r"^another_pattern_",
        ],
    },
}
```

Patterns are Python regular expressions matched against rule names.

**Customizing Classification:**

If a rule appears in the wrong stage or is not classified:
1. Check the console output for warnings about unclassified nodes
2. Add appropriate regex patterns to the relevant stage in `WORKFLOW_STAGES`
3. Re-run the script to regenerate diagrams

**CI/CD Integration:**

This script is designed to be run in CI pipelines:

```yaml
# Example: GitHub Actions
name: Update DAG Diagrams

on:
  push:
    paths:
      - 'spimquant/workflow/**'

jobs:
  update-diagrams:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'
      
      - name: Install dependencies
        run: |
          pip install snakemake snakebids
      
      - name: Generate DAG diagrams
        run: |
          python3 docs/scripts/generate_dag_diagrams.py
      
      - name: Commit updated diagrams
        run: |
          git config user.name "GitHub Actions"
          git config user.email "actions@github.com"
          git add docs/figures/*.mermaid
          git commit -m "Update DAG diagrams" || echo "No changes"
          git push
```

**Troubleshooting:**

*Issue: "No module named 'zarrnii'"*
- Install additional dependencies: `pip install zarrnii ngff-zarr dask`

*Issue: "No storage provider found"*
- Install storage plugins: `pip install snakemake-storage-plugin-http`

*Issue: Rules not classified correctly*
- Check the rule names in your Snakefile
- Update regex patterns in `WORKFLOW_STAGES`
- Rules should follow naming conventions for automatic classification

**Viewing Diagrams:**

Mermaid diagrams can be viewed directly on GitHub thanks to native Mermaid support, or using:
- [Mermaid Live Editor](https://mermaid.live/)
- Any Markdown viewer with Mermaid support
