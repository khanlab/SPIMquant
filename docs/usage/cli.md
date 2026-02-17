# Command Line Interface

The SPIMquant command-line interface is built on Snakemake and follows the BIDS App specification.

## Basic Usage

```bash
pixi run spimquant <bids_dir> <output_dir> <analysis_level> [options]
```

### Required Arguments

- `bids_dir`: Path to BIDS dataset (can be local path or S3/GCS URL)
- `output_dir`: Path to output directory
- `analysis_level`: Either `participant` or `group`

## Analysis Levels

### Participant Level

Process individual subjects:

```bash
pixi run spimquant /path/to/bids /path/to/output participant --cores all
```

This performs:

1. Template registration
2. Segmentation
3. Atlas-based quantification
4. Quality control outputs

### Group Level

Perform group statistical analysis:

```bash
pixi run spimquant /path/to/bids /path/to/output group \
  --contrast_column treatment \
  --contrast_values control drug \
  --cores all
```

Requires completed participant-level analysis.

## Common Options

### BIDS Filtering

Filter specific subjects or sessions:

```bash
# Process specific subjects
pixi run spimquant ... participant --filter_subjects 01 02 03

# Filter by SPIM data type
pixi run spimquant ... participant --filter-spim extension='ome.zarr.zip'
```

### Parallelization

Control computational resources:

```bash
# Use all available cores
pixi run spimquant ... --cores all

# Use specific number of cores
pixi run spimquant ... --cores 8

# Limit concurrent jobs
pixi run spimquant ... --jobs 4
```

### Workflow Control

```bash
# Dry run (don't execute, just plan)
pixi run spimquant ... -n

# Force re-run all steps
pixi run spimquant ... --forceall

# Re-run specific rule
pixi run spimquant ... --forcerun register_to_template

# Run until specific rule
pixi run spimquant ... --until convert_to_nifti
```

### Output Options

```bash
# Generate HTML report
pixi run spimquant ... --report

# Keep temporary files
pixi run spimquant ... --notemp

# Quiet mode
pixi run spimquant ... --quiet
```

## Advanced Options

### Template Selection

<!-- TODO: Document available templates and how to use custom templates -->

```bash
# Use specific template
pixi run spimquant ... --template gubra

# Available templates: ABAv3, gubra, MBMv3, turone, MouseIn
```

### Registration Options

<!-- TODO: Document registration parameters -->

```bash
# Specify registration stain
pixi run spimquant ... --registration_stain YOPRO
```

### Segmentation Options

<!-- TODO: Document segmentation methods and parameters -->

```bash
# Choose segmentation method
pixi run spimquant ... --segmentation_method threshold
```

### Cloud Storage

```bash
# Read from S3
pixi run spimquant s3://bucket/bids /local/output participant --cores all

# Read from GCS
pixi run spimquant gs://bucket/bids /local/output participant --cores all
```

## Group Analysis Options

When using `analysis_level group`:

```bash
pixi run spimquant /bids /output group \
  --contrast_column treatment \        # Column in participants.tsv
  --contrast_values control drug \     # Values to compare
  --cores all
```

Generates:

- `*_groupstats.tsv`: Statistical test results
- `*_groupstats.png`: Heatmap visualizations
- `*_groupstats.nii`: 3D volumetric maps
- `*_groupavgsegstats.tsv`: Group-averaged statistics
- `*_groupavg.nii.gz`: Group-averaged maps

## Configuration File

<!-- TODO: Link to configuration documentation -->

SPIMquant can also be configured via YAML files. See [Configuration Guide](configuration.md) for details.

## Environment Variables

<!-- TODO: Document relevant environment variables -->

```bash
# Set Dask configuration
export DASK_CONFIG=/path/to/dask.yaml

# AWS credentials for S3
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
```

## Snakemake Options

SPIMquant inherits all Snakemake CLI options. Key options include:

<!-- TODO: Add comprehensive Snakemake options relevant to SPIMquant -->

```bash
# Visualization
--dag                 # Output DAG visualization
--rulegraph          # Output rule graph

# Debugging
--debug              # Enable debug mode
--printshellcmds     # Print shell commands

# Resource management  
--resources mem_mb=64000  # Set memory limit
--latency-wait 60    # Wait time for file system
```

## Examples

### Basic Processing

```bash
# Process all subjects
pixi run spimquant ./bids ./output participant --cores all
```

### Advanced Workflow

```bash
# Process with custom settings
pixi run spimquant ./bids ./output participant \
  --filter_subjects 01 02 \
  --template gubra \
  --cores 16 \
  --report
```

### Cluster Execution

```bash
# Submit to SLURM cluster
pixi run spimquant ./bids ./output participant \
  --profile slurm \
  --jobs 100
```

## Getting Help

```bash
# Show help message
pixi run spimquant --help

# Show Snakemake help
pixi run spimquant --help-snakemake
```

For more detailed information, see:

- [Configuration Guide](configuration.md)
- [Workflow Guide](workflows.md)
- [Examples](../examples/workflows.md)