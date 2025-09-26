# GitHub Copilot Instructions for SPIMquant

## Project Overview

SPIMquant is a Snakemake-based BIDS app for quantitative analysis of SPIM (lightsheet) brain data. It performs automated nonlinear template registration and quantification of pathology from SPIM microscopy datasets.

### Key Features
- Deformable registration to a template
- Atlas-based quantification of pathology
- Coarse-grained and fine-grained parallelization using Snakemake and Dask
- Support for reading BIDS datasets directly from cloud-based object storage
- Support for simple and scalable cloud-based processing with Coiled

## Architecture & Framework

### Snakemake-based BIDS App
This is a Snakemake-based workflow packaged as a BIDS app. Key architectural points:

- **Main entry point**: `spimquant/run.py` using snakebids framework
- **Workflow definition**: `spimquant/workflow/Snakefile`
- **Configuration**: `spimquant/config/snakebids.yml`
- **Rules**: Organized in `spimquant/workflow/rules/` directory
- **Scripts**: Located in `spimquant/workflow/scripts/` directory

### CLI Usage
The CLI is configured via the config YAML file (`spimquant/config/snakebids.yml`). The app can be run using:

```bash
# Dry run workflow tests (recommended first step)
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant -n

# Full run with all cores
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant --cores all

# Using specific filtering for zarr files
spimquant /path/to/bids/dir /path/to/output/dir participant --filter-spim extension='ome.zarr.zip'

# Generate HTML report after processing
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant --report
```

## Development Guidelines

### Snakemake Script Patterns
When writing Snakemake Python scripts in `spimquant/workflow/scripts/`:

1. **Use the snakemake object**: Scripts automatically have access to the `snakemake` object
2. **No argument parsing needed**: Access inputs, outputs, and parameters via snakemake object
3. **Common patterns**:
   ```python
   # Access inputs
   input_file = snakemake.input.spim
   
   # Access outputs  
   output_file = snakemake.output.nii
   
   # Access parameters
   level = snakemake.params.level
   
   # Access wildcards
   stain = snakemake.wildcards.stain
   
   # Access threads
   dask.config.set(scheduler="threads", num_workers=snakemake.threads)
   ```

### BIDS Function Usage
Use the `bids()` function from snakebids for defining file paths in workflows:

```python
# In Snakemake rules
output:
    nii=bids(
        root=root,
        datatype="micr",
        stain="{stain}",
        level="{level}",
        suffix="SPIM.nii",
        **inputs["spim"].wildcards,
    ),
```

This ensures BIDS-compliant file naming and organization.

### Configuration System
- Main config: `spimquant/config/snakebids.yml`
- Defines BIDS inputs, analysis levels, and app-specific parameters
- Example config templates in `examples/` directory
- Users copy and modify config files for their specific datasets

### Testing & Quality Assurance

#### Dry Run Testing
Always test workflows with the `-n` flag before full execution:
```bash
spimquant /path/to/bids/dir /path/to/output/dir participant -n
```

#### Quality Checks
The project uses automated formatting and linting:
```bash
pixi run quality_check  # Check code formatting
pixi run quality_fix    # Fix code formatting
```

Tools used:
- **isort**: Python import sorting
- **black**: Python code formatting  
- **snakefmt**: Snakemake file formatting

### Project Structure
```
spimquant/
├── run.py                    # Main CLI entry point
├── config/
│   └── snakebids.yml        # Main configuration
├── workflow/
│   ├── Snakefile           # Main workflow definition
│   ├── rules/              # Modular rule definitions
│   └── scripts/            # Python scripts for rules
├── resources/              # Reference data and templates
└── profiles/              # Execution profiles
```

## Documentation Requirements

### For New Features
When adding new features, ensure:

1. **Update documentation**: Add usage examples and parameter descriptions
2. **Update config schema**: Document new configuration options in `snakebids.yml`
3. **Add tests**: Include dry-run tests for new workflows
4. **Code comments**: Document complex algorithms and data processing steps
5. **README updates**: Add new features to the features list and usage examples

### Code Documentation Standards
- Use docstrings for functions and classes
- Comment complex Snakemake rules with purpose and parameter explanations
- Document file format expectations and outputs
- Include references for scientific methods and algorithms

## Development Environment

### Setup
```bash
# Install pixi (dependency manager)
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and setup
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant
pixi install --environment dev  # For development with linters
```

### Hardware Requirements
- Sufficient memory (at least 16GB+, preferably more)
- Multiple cores benefit parallelization
- Fast local storage for temporary files
- For `greedy` registration: significant memory during template registration

### Dependencies
Managed through pixi/conda:
- **Core**: snakemake, snakebids, dask, zarr
- **Image processing**: antspyx, greedyreg, zarrnii
- **Cloud storage**: s3fs, gcsfs, storage plugins
- **Development**: black, isort, snakefmt, pytest

## Specific Implementation Notes

### SPIM Data Handling
- Primary data format: OME-Zarr files
- Support for cloud storage (S3, GCS)
- Multi-level/multi-resolution data processing
- Channel-specific processing for different stains

### Registration Workflow
- Uses `greedy` for deformable registration
- Template-based registration with configurable templates:
  - ABAv3 (Allen Brain Atlas v3)
  - gubra
  - MBMv3 (Mouse Brain Maps v3)
  - turone
  - MouseIn (for MRI registration)
- Atlas-based segmentation and quantification
- Memory-intensive operations require careful resource management
- Supports multi-level/multi-resolution processing for efficiency

### Segmentation & Quantification
- Configurable stains for registration: PI, YOPRO, AutoF, etc.
- Segmentation stains: abeta, Abeta, BetaAmyloid, AlphaSynuclein, Iba1, ChAT
- Segmentation methods:
  - `threshold`: Simple threshold-based segmentation
  - `otsu+k3i2`: Otsu thresholding with k-means clustering
- Intensity correction methods: gaussian, n4
- Field fraction calculations for quantitative analysis

### Parallelization Strategy
- Snakemake handles workflow-level parallelization
- Dask for array processing within scripts
- SLURM integration for cluster execution
- Thread configuration through snakemake.threads

### Execution Profiles
- Default profile available in `spimquant/profiles/default/`
- Configurable execution environments for different compute resources
- Cloud execution support through Coiled integration
- SLURM executor plugin for cluster environments

## Best Practices

1. **Always dry-run first**: Use `-n` flag to validate workflows
2. **Use BIDS functions**: Consistent file naming and organization
3. **Leverage snakemake object**: No need for custom argument parsing
4. **Document new features**: Code, configuration, and user documentation
5. **Test incrementally**: Small changes with dry-run validation
6. **Follow code style**: Use provided linting tools
7. **Memory awareness**: Monitor resource usage for large datasets
8. **Cloud-ready**: Consider storage plugins for remote data access