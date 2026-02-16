# Frequently Asked Questions

<!-- TODO: Review this section for accuracy - some content needs revision and verification -->

## General Questions

### What is SPIMquant?

SPIMquant is a Snakemake-based BIDS app for quantitative analysis of SPIM (lightsheet) brain microscopy data. It performs automated registration to brain templates and atlas-based quantification of pathology.

### What data formats does SPIMquant support?

SPIMquant works with BIDS-formatted datasets containing OME-Zarr files. We recommend using [SPIMprep](https://github.com/khanlab/SPIMprep) to convert your raw microscopy data to BIDS format.

### What templates are supported?

SPIMquant supports:
- ABAv3 (Allen Brain Atlas v3)
- gubra
- MBMv3 (Marmoset Brain Maps v3)
- turone
- MouseIn (for MRI registration)

### Can I use custom templates?

Yes! <!-- TODO: Add link to custom template documentation --> See the [Custom Templates Tutorial](tutorials/custom_templates.md) for details.

## Installation

### What operating systems are supported?

SPIMquant currently supports Linux (64-bit). Windows and macOS support may be added in future releases.

### How much memory do I need?

Minimum 32 GB RAM is recommended. The `greedy` registration step can require 16-32 GB for a single job. See [Hardware Requirements](getting_started/hardware.md) for details.

### Do I need a GPU?

No, SPIMquant does not require a GPU. All processing is CPU-based.

### How do I install dependencies?

Dependencies are managed automatically through pixi. Simply run `pixi install` after cloning the repository. See [Installation Guide](getting_started/installation.md).

## Usage

### How long does processing take?

Processing time varies significantly based on:
- Dataset size and resolution
- Available compute resources
- Number of subjects
- Specific workflow steps

Typical times: 2-4 hours per subject on a 16-core machine with 64 GB RAM.

### Can I process data in the cloud?

Yes! SPIMquant supports reading from S3 and GCS, and can be deployed on cloud infrastructure. See [Cloud Processing Guide](usage/cloud.md).

### How do I resume after a failure?

Simply re-run the same command. Snakemake automatically detects completed steps and resumes from the last successful job.

### Can I process multiple subjects in parallel?

Yes, use `--cores all` to utilize all available cores. Snakemake will automatically parallelize across subjects.

### How do I reduce memory usage?

- Use fewer cores: `--cores 4` instead of `--cores all`
- Process subjects sequentially
- Use lower resolution pyramid levels
- Reduce concurrent jobs

## Data and Formats

### What is OME-Zarr?

OME-Zarr is a cloud-optimized format for storing multidimensional microscopy data. It supports efficient storage, streaming, and parallel access.

### How do I convert my data to BIDS format?

Use [SPIMprep](https://github.com/khanlab/SPIMprep) to convert raw SPIM data to BIDS-formatted OME-Zarr files.

### Can I use compressed Zarr files?

Yes, SPIMquant supports Zarr zipstores (`.ome.zarr.zip`). Use the filter option:
```bash
--filter-spim extension='ome.zarr.zip'
```

### What stains are supported for registration?

Common registration stains include:
- PI (Propidium Iodide)
- YOPRO
- AutoF (Autofluorescence)

### What stains can be segmented?

Segmentation is supported for:
- Abeta/BetaAmyloid
- AlphaSynuclein
- Iba1
- ChAT
- Custom stains

## Registration and Segmentation

### Why did registration fail?

Common causes:
- Poor data quality
- Template mismatch
- Insufficient memory
- Incorrect stain specification

Check input data quality and ensure sufficient memory is available.

### How do I choose segmentation methods?

- **threshold**: Simple intensity thresholding, fast
- **otsu+k3i2**: Otsu + k-means, more robust

<!-- TODO: Add link to segmentation guide --> See [Segmentation Methods](howto/segmentation.md) for details.

### Can I adjust registration parameters?

Yes, registration parameters can be adjusted in the configuration file. <!-- TODO: Add link to config guide --> See [Configuration Guide](usage/configuration.md).

## Group Analysis

### How many subjects do I need per group?

Minimum 3 subjects per group recommended. Larger groups (n≥10) provide better statistical power.

### What statistical tests are used?

Currently t-tests for two-group comparisons. <!-- TODO: Update when more methods added -->

### How do I interpret effect sizes?

Cohen's d interpretation:
- Small: d < 0.2
- Medium: 0.2 ≤ d < 0.8
- Large: d ≥ 0.8

### Can I include covariates?

<!-- TODO: Update when covariate support is added --> Covariate support is planned for future releases.

## Outputs

### Where are my results?

Results are in the output directory under `spimquant/`:
```
output/spimquant/
├── sub-01/
│   └── micr/
│       ├── *_space-template_SPIM.nii.gz
│       ├── *_dseg.nii.gz
│       └── *_segstats.tsv
```

### What format are the statistics in?

Statistics are saved as TSV (tab-separated values) files that can be opened in Excel, R, Python, or other analysis tools.

### How do I visualize 3D results?

Use neuroimaging software such as:
- FSLeyes
- ITK-SNAP
- 3D Slicer
- Imaris

### Can I generate reports?

Yes, use the `--report` flag to generate an HTML report:
```bash
pixi run spimquant ... --report
```

## Performance and Optimization

### How do I speed up processing?

- Use more CPU cores
- Use SSD storage for temporary files
- Process lower resolution first for testing
- Use cloud resources for scaling

### Why is processing slow?

Common bottlenecks:
- Slow storage (HDD vs SSD)
- Limited memory causing swapping
- Network storage latency
- Large dataset size

### How much disk space do I need?

Plan for 5-10x input data size for intermediate files. Use `--delete-temp-output` to clean up after completion.

## Troubleshooting

### "Out of memory" errors

- Reduce concurrent jobs: `--cores 4`
- Close other applications
- Process subjects one at a time
- Increase available RAM

### "File not found" errors

- Check BIDS naming conventions
- Verify file extensions
- Ensure data is in correct directory structure
- Check file permissions

### "Lock directory exists" errors

Previous run was interrupted. Remove locks:
```bash
rm -rf /output/.snakemake/locks/
```

### Snakemake errors

- Run with `-n` for dry run to identify issues
- Use `--verbose` for detailed output
- Check `.snakemake/log/` for error logs
- Ensure all dependencies are installed

## Getting Help

### Where can I get help?

1. Check this FAQ
2. Search [GitHub Issues](https://github.com/khanlab/SPIMquant/issues)
3. Ask on [GitHub Discussions](https://github.com/khanlab/SPIMquant/discussions)
4. Open a new issue with detailed information

### What information should I include in bug reports?

Include:
- SPIMquant version
- Operating system and version
- Full error message and stack trace
- Command used
- Relevant log files
- Steps to reproduce

### How can I contribute?

See the [Contributing Guide](contributing.md) for how to contribute code, documentation, or report issues.

## Future Features

### What features are planned?

<!-- TODO: Update with roadmap -->

Planned features include:
- Additional statistical methods
- More segmentation algorithms
- Enhanced visualization tools
- Improved cloud integration
- Multi-modal registration

### Can I request features?

Yes! Open a [feature request](https://github.com/khanlab/SPIMquant/issues/new) on GitHub.