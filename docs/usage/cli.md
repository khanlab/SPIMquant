# Command Line Interface

SPIMquant provides a comprehensive command-line interface built on Snakemake and Snakebids frameworks. This page covers all available options and provides practical usage examples.

## Core BIDS App Interface

### Basic Usage Pattern

```bash
spimquant <bids_dir> <output_dir> <analysis_level> [OPTIONS]
```

**Required Arguments:**
- `bids_dir`: Path or URI to BIDS dataset (local path, s3://, gs://, etc.)
- `output_dir`: Local output directory path
- `analysis_level`: Processing level (currently only 'participant' supported)

### Core Arguments and Options

```{argparse}
:ref: spimquant.run.get_parser
:prog: spimquant
```

## Common Usage Examples

### Basic Processing

**Dry Run (Recommended First Step):**
```bash
# Check workflow without executing
spimquant /path/to/bids/data /path/to/output participant \
  --dry-run --cores 1
```

**Standard Processing:**
```bash
# Process all subjects using containers
spimquant /path/to/bids/data /path/to/output participant \
  --cores all --use-apptainer
```

**Conda Environment:**
```bash
# Use conda for dependencies instead of containers
spimquant /path/to/bids/data /path/to/output participant \
  --cores 8 --use-conda
```

### Template and Atlas Selection

**Choose Template:**
```bash
# Use Gubra template instead of default ABAv3
spimquant /path/to/data /path/to/output participant \
  --template gubra --cores 8
```

**Select Atlas Resolution:**
```bash
# Use specific atlas parcellation
spimquant /path/to/data /path/to/output participant \
  --atlas_segs roi82 --cores 8

# Use multiple atlas levels
spimquant /path/to/data /path/to/output participant \
  --atlas_segs roi22 roi82 roi198 --cores 8
```

**Hemisphere-Specific Analysis:**
```bash
# Process only left hemisphere
spimquant /path/to/data /path/to/output participant \
  --template_crop left --cores 8
```

### Stain Configuration

**Custom Registration Stains:**
```bash
# Specify stains for template registration
spimquant /path/to/data /path/to/output participant \
  --stains_for_reg DAPI PI AutoF --cores 8
```

**Custom Analysis Stains:**
```bash
# Specify stains for quantification
spimquant /path/to/data /path/to/output participant \
  --stains_for_seg abeta pTau Iba1 GFAP --cores 8
```

### Processing Parameters

**Resolution Control:**
```bash
# High-resolution processing
spimquant /path/to/data /path/to/output participant \
  --registration_level 3 --segmentation_level 0 --cores 4

# Fast/low-resolution processing  
spimquant /path/to/data /path/to/output participant \
  --registration_level 6 --segmentation_level 3 --cores 8 --sloppy
```

**Segmentation Method:**
```bash
# Use manual threshold
spimquant /path/to/data /path/to/output participant \
  --seg_method threshold --seg_threshold 85 --cores 8

# Use adaptive Otsu thresholding
spimquant /path/to/data /path/to/output participant \
  --seg_method otsu+k3i2 --cores 8
```

### Cloud Data Processing

**Amazon S3:**
```bash
# Read from S3 bucket
spimquant s3://my-bucket/bids-dataset /local/output participant \
  --cores 16 --use-apptainer

# With custom credentials
export AWS_PROFILE=my-profile
spimquant s3://my-bucket/bids-dataset /local/output participant \
  --cores 16
```

**Google Cloud Storage:**
```bash
# Read from GCS bucket
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json
spimquant gs://my-bucket/bids-dataset /local/output participant \
  --cores 16 --use-apptainer
```

### Subject Selection

**Process Specific Subjects:**
```bash
# Single subject
spimquant /path/to/data /path/to/output participant \
  --participant_label 001 --cores 8

# Multiple subjects
spimquant /path/to/data /path/to/output participant \
  --participant_label 001 002 005 --cores 8

# Exclude subjects
spimquant /path/to/data /path/to/output participant \
  --exclude_participant_label 003 004 --cores 8
```

### Data Format Options

**OME-Zarr Zip Files:**
```bash
# For .ome.zarr.zip archives
spimquant /path/to/data /path/to/output participant \
  --filter-spim extension='ome.zarr.zip' --cores 8
```

**Custom BIDS Filters:**
```bash
# Filter by acquisition parameter
spimquant /path/to/data /path/to/output participant \
  --filter-spim acquisition='highres' --cores 8

# Multiple filters
spimquant /path/to/data /path/to/output participant \
  --filter-spim sample='brain' staining='immunofluor' --cores 8
```

## Resource Management

### CPU and Memory Control

**Core Allocation:**
```bash
# Use all available cores
spimquant /path/to/data /path/to/output participant --cores all

# Limit to 8 cores
spimquant /path/to/data /path/to/output participant --cores 8

# Single-threaded processing
spimquant /path/to/data /path/to/output participant --cores 1
```

**Memory Limits:**
```bash
# Limit memory per job
spimquant /path/to/data /path/to/output participant \
  --cores 4 --resources mem_mb=8000

# Total workflow memory limit
spimquant /path/to/data /path/to/output participant \
  --cores 8 --resources mem_mb=16000 total_mem=32000
```

### Cluster Execution

**SLURM Cluster:**
```bash
# Submit to SLURM
spimquant /path/to/data /path/to/output participant \
  --executor slurm --jobs 10 \
  --slurm-partition compute \
  --slurm-account my_account
```

**Custom Cluster Profile:**
```bash
# Use cluster profile
spimquant /path/to/data /path/to/output participant \
  --profile cluster_profile --jobs 20
```

## Advanced Options

### Workflow Control

**Specific Target Rules:**
```bash
# Run only registration
spimemake /path/to/data /path/to/output participant \
  --target all_templatereg --cores 8

# Run only segmentation
spimquant /path/to/data /path/to/output participant \
  --target all_segment --cores 8
```

**Workflow Continuation:**
```bash
# Continue from where left off
spimquant /path/to/data /path/to/output participant \
  --cores 8 --rerun-incomplete

# Force re-run specific steps
spimquant /path/to/data /path/to/output participant \
  --cores 8 --forcerun templatereg
```

### Debugging and Development

**Verbose Output:**
```bash
# Increase verbosity
spimquant /path/to/data /path/to/output participant \
  --cores 8 --verbose

# Debug mode
spimquant /path/to/data /path/to/output participant \
  --cores 8 --debug-dag
```

**Keep Intermediate Files:**
```bash
# Don't delete temporary files
spimquant /path/to/data /path/to/output participant \
  --cores 8 --notemp
```

## Snakemake Integration

SPIMquant inherits all Snakemake command-line options. Any Snakemake option can be used alongside SPIMquant-specific options.

### Core Snakemake Options

```{argparse}
:module: snakebids.snakemake_compat
:func: get_argument_parser
:prog: snakemake
```

### Useful Snakemake Options

**Workflow Inspection:**
```bash
# Show workflow DAG
spimquant /path/to/data /path/to/output participant \
  --dag | dot -Tpng > workflow.png

# List all rules
spimquant /path/to/data /path/to/output participant \
  --list-rules

# Show rule dependencies
spimquant /path/to/data /path/to/output participant \
  --rulegraph | dot -Tpng > rules.png
```

**Reporting and Monitoring:**
```bash
# Generate HTML report
spimquant /path/to/data /path/to/output participant \
  --cores 8 --report report.html

# Real-time monitoring
spimquant /path/to/data /path/to/output participant \
  --cores 8 --printshellcmds
```

## Configuration Files

### Using Custom Configuration

```bash
# Specify custom config file
spimquant /path/to/data /path/to/output participant \
  --configfile custom_config.yml --cores 8

# Override specific config values
spimquant /path/to/data /path/to/output participant \
  --config template=gubra registration_level=4 --cores 8
```

### Example Configuration Override

```bash
# Multiple configuration overrides
spimquant /path/to/data /path/to/output participant \
  --config \
    template=ABAv3 \
    registration_level=5 \
    segmentation_level=2 \
    correction_method=n4 \
    seg_method=otsu+k3i2 \
  --cores 8
```

## Common Command Combinations

### Production Processing

```bash
# Robust production run
spimquant /path/to/data /path/to/output participant \
  --use-apptainer \
  --cores all \
  --resources mem_mb=8000 \
  --template ABAv3 \
  --atlas_segs roi82 roi198 \
  --correction_method n4 \
  --seg_method otsu+k3i2 \
  --report production_report.html
```

### Development Testing

```bash
# Quick test run
spimquant /path/to/data /path/to/output participant \
  --participant_label 001 \
  --sloppy \
  --registration_level 6 \
  --segmentation_level 4 \
  --cores 4 \
  --dry-run
```

### High-Quality Analysis

```bash
# Maximum quality processing
spimquant /path/to/data /path/to/output participant \
  --use-apptainer \
  --cores 16 \
  --registration_level 3 \
  --segmentation_level 0 \
  --correction_method n4 \
  --template ABAv3 \
  --atlas_segs all roi22 roi82 roi198 \
  --resources mem_mb=16000
```

