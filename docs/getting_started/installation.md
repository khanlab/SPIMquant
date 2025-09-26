## SPIMquant

SPIMquant is a Snakebids BIDS App for quantitative analysis of SPIM (lightsheet) brain microscopy data. It performs automated nonlinear template registration and quantification of pathology from SPIM microscopy datasets.

### Key Features

- **Deformable registration** to anatomical templates (ABAv3, Gubra, MBMv3, Turone)
- **Atlas-based quantification** of pathology and microscopic features
- **Multi-stain support** for different fluorescent markers
- **Scalable processing** with Snakemake and containerized workflows
- **Cloud support** for reading BIDS datasets from object storage (S3, GCS)
- **Flexible configuration** for different experimental setups

### Hardware Requirements

- **Memory**: At least 16GB RAM, preferably 32GB+ for full-resolution processing
  - The deformable registration (`greedy`) can be memory-intensive
  - Memory usage scales with image resolution and number of parallel jobs
- **Storage**: 
  - Fast local storage recommended for temporary files (SSD preferred)
  - Sufficient disk space for intermediate files (can be 5-10x input data size)
- **CPU**: Multi-core processor recommended for parallel processing

### Software Requirements

**Recommended Setup (Linux with containers):**
- Linux operating system (Ubuntu 20.04+ or similar)
- Python 3.11+
- Singularity/Apptainer for containerized dependencies
- Git for installation

**Alternative Setup (native dependencies):**
- `itk-snap` (includes `c3d` command-line tools)
- `greedy` registration toolkit
- `ANTs` for additional image processing
- Python environment with dependencies from `pyproject.toml`

### Sample Dataset

A sample dataset is available for testing and learning:

[Download Sample Dataset](https://drive.google.com/file/d/1-eVG_1VREKCE8auI81jyW4onyUcpzXI7/view?usp=sharing)

The sample dataset contains:
- OME-Zarr formatted SPIM microscopy data in BIDS structure
- Downsampled mouse brain scans with multiple stains
- Reference template files for registration
- Example configuration files

## Installation

### Method 1: Direct Installation from GitHub (Recommended)

Install SPIMquant directly from the GitHub repository:

```bash
pip install git+https://github.com/khanlab/spimquant.git
```

This will install the latest version with all Python dependencies.

### Method 2: Development Installation

If you plan to modify the code or contribute to development:

```bash
# Clone the repository
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant

# Install in development mode
pip install -e .
```

### Method 3: Using Pixi (For Developers)

SPIMquant includes a `pixi.toml` configuration for reproducible environments:

```bash
# Install pixi if not already installed
curl -fsSL https://pixi.sh/install.sh | bash

# Clone and setup environment
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant

# Create environment and install
pixi install
pixi shell  # Activate environment
```

### Container Dependencies

SPIMquant uses containerized tools for image processing. These containers are automatically downloaded when needed:

- **c3d**: Image conversion and basic operations
- **greedy**: Deformable image registration  
- **ANTs**: Advanced image processing and registration

### Verification

Test your installation:

```bash
# Check if spimquant is available
spimquant --help

# Verify version
spimquant --version
```

### Configuration Setup

After installation, you need to prepare a configuration file:

1. **Copy the template configuration:**
   ```bash
   mkdir -p ~/.config/spimquant
   wget https://github.com/khanlab/SPIMquant/raw/main/spimquant/config/snakebids.yml \
        -O ~/.config/spimquant/snakebids.yml
   ```

2. **Or use the built-in default configuration** (recommended for beginners)

The configuration file controls processing parameters, template choices, and analysis options. See the [Configuration section](../usage/configuration.md) for detailed parameter descriptions.

## Basic Usage

SPIMquant follows the BIDS App specification, requiring a BIDS-formatted dataset as input. The [SPIMprep](https://github.com/khanlab/SPIMprep) workflow can convert raw microscopy data into BIDS format.

### Quick Start

1. **Dry run** (recommended first step to check configuration):
   ```bash
   spimquant /path/to/bids/dataset /path/to/output participant --dry-run --cores 1
   ```

2. **Full processing** using all available CPU cores:
   ```bash
   spimquant /path/to/bids/dataset /path/to/output participant --cores all --use-conda
   ```

3. **Processing with containers** (recommended for reproducibility):
   ```bash
   spimquant /path/to/bids/dataset /path/to/output participant --cores all --use-apptainer
   ```

### Input Requirements

**BIDS Dataset Structure:**
```
dataset/
├── sub-001/
│   └── micr/
│       ├── sub-001_sample-brain_stain-PI_SPIM.ome.zarr/
│       ├── sub-001_sample-brain_stain-abeta_SPIM.ome.zarr/
│       └── sub-001_sample-brain_stain-Iba1_SPIM.ome.zarr/
├── sub-002/
│   └── micr/
│       └── ...
└── dataset_description.json
```

**Supported File Formats:**
- **OME-Zarr**: `*.ome.zarr` (directories) or `*.ome.zarr.zip` (archives)
- **Multiple stains**: Different fluorescent markers per file
- **Multi-resolution**: Pyramid structures supported for efficient processing

### Special Options

**For Zarr Zip Archives:**
```bash
spimquant /path/to/bids/dataset /path/to/output participant \
  --filter-spim extension='ome.zarr.zip' \
  --cores all --use-apptainer
```

**Cloud Data Processing:**
```bash
# Read from S3 bucket
spimquant s3://bucket-name/bids-dataset /path/to/output participant --cores all

# Read from Google Cloud Storage
spimquant gs://bucket-name/bids-dataset /path/to/output participant --cores all
```

**Memory-Conscious Processing:**
```bash
# Reduce memory usage for large datasets
spimquant /path/to/bids/dataset /path/to/output participant \
  --registration_level 6 \
  --segmentation_level 2 \
  --cores 4 --resources mem_mb=8000
```


## Output and Reporting

### Processing Output

SPIMquant generates a comprehensive set of outputs in BIDS-derivatives format:

```
output/
├── spimquant/
│   ├── sub-001/
│   │   └── micr/
│   │       ├── sub-001_space-ABAv3_stain-PI_SPIM.nii         # Registered images
│   │       ├── sub-001_space-ABAv3_seg-all_fieldfrac.nii    # Field fraction maps
│   │       ├── sub-001_seg-all_stain-abeta_blobdensity.nii  # Quantification results
│   │       └── ...
│   ├── tpl-ABAv3/                                           # Template space files
│   │   ├── tpl-ABAv3_desc-LR_dseg.nii.gz                   # Atlas segmentation
│   │   └── tpl-ABAv3_desc-LR_dseg.tsv                      # Region labels
│   └── logs/                                                # Processing logs
└── work/                                                    # Temporary files
```

### Generating HTML Reports

Create a detailed processing report with Snakemake's built-in reporting:

```bash
spimquant /path/to/bids/dataset /path/to/output participant --report report.html
```

The report includes:
- **Processing graph**: Visual workflow representation
- **Job details**: Runtime, memory usage, and commands
- **File provenance**: Input/output relationships
- **Configuration**: All parameters used
- **Performance metrics**: Timing and resource usage

### Quality Control

Inspect results using image viewers:

```bash
# View registered images
itksnap output/spimquant/sub-001/micr/sub-001_space-ABAv3_stain-PI_SPIM.nii

# Load template and segmentation
itksnap \
  -g output/spimquant/tpl-ABAv3/tpl-ABAv3_anat.nii.gz \
  -s output/spimquant/tpl-ABAv3/tpl-ABAv3_desc-LR_dseg.nii.gz
```


