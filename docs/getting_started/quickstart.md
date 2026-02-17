# Quick Start

This guide will walk you through running your first SPIMquant workflow.

## Prerequisites

Before starting, ensure you have:

1. [Installed SPIMquant](installation.md)
2. A BIDS-formatted SPIM dataset (see [SPIMprep](https://github.com/khanlab/SPIMprep) for converting raw data)
3. Sufficient disk space and memory (see [Hardware Requirements](hardware.md))

## Step 1: Prepare Your Data

SPIMquant requires a BIDS dataset with SPIM microscopy data. Your directory structure should look like:

```
bids_dataset/
├── dataset_description.json
├── participants.tsv
└── sub-01/
    └── micr/
        └── sub-01_sample-brain_stain-YOPRO_SPIM.ome.zarr.zip
```

!!! tip "Using SPIMprep"
    Use [SPIMprep](https://github.com/khanlab/SPIMprep) to convert your raw microscopy data to BIDS format.

## Step 2: Run a Dry Run

Always start with a dry run to verify the workflow:

```bash
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant -n
```

The `-n` flag performs a dry run that:

- Validates your input data
- Shows what jobs would be executed
- Checks for missing files or configuration issues
- Does not actually process any data

!!! tip "Using Test Data"
    For your first run, use the absolute path to the test dataset:
    ```bash
    pixi run spimquant /absolute/path/to/SPIMquant/tests/bids_ds /tmp/output participant -n
    ```

## Step 3: Run the Workflow

Once the dry run succeeds, run the full workflow:

```bash
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant --cores all
```

Options explained:

- `participant`: Run subject-level analysis
- `--cores all`: Use all available CPU cores
- Alternative: `--cores 8` to use 8 cores

### For Zarr Zipstore Data

If your SPIM files end with `*.ome.zarr.zip`, add the filter option:

```bash
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant \
  --filter-spim extension='ome.zarr.zip' \
  --cores all
```

## Step 4: Monitor Progress

SPIMquant will display progress as it processes:

```
Building DAG of jobs...
Job counts:
    count   jobs
    1       all
    3       create_template_crop
    ...
[Mon Jan 15 10:23:45 2024]
rule convert_to_nifti:
    input: ...
    output: ...
```

Processing time depends on:

- Dataset size
- Available compute resources
- Number of subjects
- Resolution of data

## Step 5: View Results

After processing completes, find your results in the output directory:

```
output/
└── spimquant/
    ├── sub-01/
    │   └── micr/
    │       ├── sub-01_space-template_SPIM.nii.gz
    │       ├── sub-01_space-template_dseg.nii.gz
    │       └── sub-01_segstats.tsv
    └── qc/
        └── sub-01_registration_overlay.png
```

Key output files:

- `*_space-template_SPIM.nii.gz`: Registered SPIM data
- `*_dseg.nii.gz`: Segmentation results
- `*_segstats.tsv`: Quantitative statistics by brain region

## Step 6: Generate a Report

Create an HTML report of the workflow:

```bash
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant --report
```

This generates `report.html` with:

- Workflow graph with all jobs
- Runtime statistics
- Configuration details
- Clickable nodes to inspect code

## Common Issues

### Memory Errors

If you encounter out-of-memory errors:

```bash
# Use fewer cores
pixi run spimquant ... --cores 4

# Or process one subject at a time
pixi run spimquant ... --cores all --until convert_to_nifti
```

### Input Data Not Found

Ensure your data follows BIDS naming conventions:

- Files must have proper BIDS entities (sub-, sample-, stain-)
- Must be in correct directory structure (sub-XX/micr/)
- Check file extensions match your filter

### Registration Failures

If registration fails:

1. Check input data quality
2. Verify template compatibility
3. Try different registration stain
4. Adjust registration parameters in config

## Next Steps

- [Configuration Guide](../usage/configuration.md): Customize your workflow
- [Group Analysis](../usage/group_analysis.md): Compare across experimental groups  
- [CLI Reference](../usage/cli.md): Explore all command options
- [Tutorials](../tutorials/basic_registration.md): In-depth workflow examples

## Example Workflows

### Basic Processing

```bash
# Process all subjects with default settings
pixi run spimquant ./bids ./output participant --cores all
```

### Custom Template

```bash
# Use specific template
pixi run spimquant ./bids ./output participant \
  --template gubra \
  --cores all
```

### Cloud Data

```bash
# Process data stored in S3
pixi run spimquant s3://bucket/bids ./output participant \
  --cores all
```

### Group Statistics

```bash
# After participant-level, run group analysis
pixi run spimquant ./bids ./output group \
  --contrast_column treatment \
  --contrast_values control drug \
  --cores all
```