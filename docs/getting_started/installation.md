# Installation

This guide will help you install SPIMquant and its dependencies.

## Prerequisites

### Hardware Requirements

SPIMquant is computationally intensive and requires:

- **Memory**: At least 16GB RAM (32GB+ recommended)
  - The `greedy` deformable registration can consume significant memory
  - More memory allows for faster parallel processing
- **Storage**: Fast local storage for temporary files
  - Processing generates large intermediate files
  - SSD storage recommended for performance
- **CPU**: Multiple cores benefit parallelization
  - More cores = faster processing
  - Recommended: 8+ cores

### Software Requirements

- **Operating System**: Linux (64-bit)
- **pixi**: Package manager for managing dependencies
  - Manages Python packages and external tools automatically
  - Handles tools like `greedy`, `ANTs`, `c3d`

## Installation Steps

There are two ways to install SPIMquant: as a conda package, if you only
want to run it, or from a clone of the repository, if you want to modify it
or use the development tools.

### Option A: Install the conda package

Install [pixi](https://pixi.prefix.dev/latest/installation/) (see step 1
below), then:

```bash
pixi global install \
  --channel https://prefix.dev/conda-forge \
  --channel https://prefix.dev/bioconda \
  --channel https://prefix.dev/akhanf/khanlab \
  spimquant

spimquant --help
```

The package is also installable with conda, mamba or micromamba from the
same three channels.

### Option B: Install from the repository

### 1. Install pixi

If you don't have pixi installed, install it using:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

After installation, restart your terminal or source your shell configuration:

```bash
source ~/.bashrc  # or ~/.zshrc for zsh
```

### 2. Clone the Repository

```bash
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant
```

### 3. Install Dependencies

For regular usage:

```bash
pixi install
```

For development (includes formatters, linters, and visualization tools):

```bash
pixi install --environment dev
```

This will create a conda environment and install all required dependencies including:

- Python packages (snakemake, snakebids, zarr, dask, etc.)
- Image processing tools (greedy, ANTs, c3d)
- Cloud storage plugins (for S3, GCS)

## Verify Installation

Test that SPIMquant is installed correctly:

```bash
pixi run spimquant --help
```

You should see the command-line interface help message.

## Sample Dataset

<!-- TODO: Update this link when sample dataset is finalized -->

A sample dataset for testing is available [here](https://drive.google.com/file/d/1-eVG_1VREKCE8auI81jyW4onyUcpzXI7/view?usp=sharing).

The sample includes:

- Down-sampled mouse brain OME-Zarr data in BIDS format
- Reference template files
- Example configuration file

## Next Steps

- [Quick Start Guide](quickstart.md): Run your first SPIMquant workflow
- [Configuration](../usage/configuration.md): Learn about configuration options
- [CLI Reference](../usage/cli.md): Explore command-line options

## Troubleshooting

### Common Issues

**Installation fails with conda/mamba errors:**

- Try updating pixi: `curl -fsSL https://pixi.sh/install.sh | bash`
- Clear pixi cache: `rm -rf ~/.pixi/cache`

**"Command not found" after installation:**

- Ensure pixi is in your PATH
- Restart your terminal
- Use `pixi run spimquant` instead of just `spimquant`

**Out of memory errors:**

- Reduce the number of parallel jobs: `--cores 4` instead of `--cores all`
- Process smaller regions or lower resolution data first
- Increase available system memory

### Getting Help

If you encounter issues:

1. Check the [FAQ](../faq.md)
2. Search [GitHub Issues](https://github.com/khanlab/SPIMquant/issues)
3. Open a new issue with:
   - Your system information (`uname -a`)
   - SPIMquant version (`pixi run spimquant --version`)
   - Full error message and stack trace
