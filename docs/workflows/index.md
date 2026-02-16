# Workflow Documentation

SPIMquant is organized into modular Snakemake workflows, each handling a specific aspect of the SPIM data processing pipeline. The workflows are designed to be composable and can be executed independently or as part of the complete pipeline.

## Workflow Architecture

The SPIMquant workflow is divided into several specialized modules:

```
SPIM Data (OME-Zarr)
    ↓
[Import & Format Conversion] → NIfTI files at various resolutions
    ↓
[Brain Masking] → Brain tissue masks
    ↓
[Template Registration] → Spatial normalization to atlas
    ↓
[Segmentation & Quantification] → Detect pathology & compute metrics
    ↓
[Group Statistics] → Population-level analysis
```

## Core Workflow Modules

### [Import and Format Conversion](import.md)
Handles data import and format conversion between OME-Zarr and NIfTI formats.
- Template and atlas import
- Multi-resolution pyramid handling
- BIDS organization

### [Brain Masking](masking.md)
Creates brain tissue masks using Gaussian mixture modeling with template priors.
- Tissue classification (Atropos)
- Template-guided masking
- Background removal

### [Template Registration](templatereg.md)
Registers subject SPIM data to standard template space.
- Intensity correction (N4)
- Affine and deformable registration
- Forward/inverse transformations
- Quality control reporting

### [Segmentation and Quantification](segmentation.md)
Segments pathology markers and computes quantitative metrics.
- Bias field correction
- Multi-channel segmentation
- Region properties extraction
- Atlas-based quantification
- Colocalization analysis

### [Group Statistics](groupstats.md)
Performs group-level statistical analysis across subjects.
- Group comparisons
- Statistical testing per brain region
- Heatmap visualization
- Volumetric statistical maps

## Optional Workflow Modules

### [Patch Extraction](patches.md)
Extracts 3D patches for visualization and machine learning.
- Atlas-guided sampling
- Multiple data types (raw, corrected, masks)
- Imaris dataset export

### [MRI Preprocessing](preproc_mri.md)
Co-registers in-vivo MRI with ex-vivo SPIM data (when available).
- MRI preprocessing
- Cross-modality registration
- Tissue deformation analysis

### [Common Utilities](common.md)
Shared utility functions used across all workflows.
- Path management
- BIDS helpers
- Configuration accessors

## Data Flow and Dependencies

The workflows are organized in a dependency hierarchy:

1. **Import** must run first to provide reference data
2. **Masking** requires imported templates and subject data
3. **Registration** requires masks and templates
4. **Segmentation** requires registration outputs for atlas-based analysis
5. **Group Statistics** requires multiple subjects' segmentation results

Snakemake automatically resolves these dependencies and executes rules in the correct order.

## File Formats

### Inputs
- **OME-Zarr** (`.ome.zarr`): Multi-resolution SPIM microscopy data
- **NIfTI** (`.nii.gz`): Template anatomies and atlases
- **TSV** (`.tsv`): Atlas label tables and metadata

### Outputs
- **NIfTI** (`.nii.gz`): Processed images, masks, statistical maps
- **Parquet** (`.parquet`): Region properties tables (efficient for large datasets)
- **TSV** (`.tsv`): Statistical summaries, per-region metrics
- **PNG** (`.png`): Quality control visualizations
- **HTML** (`.html`): Interactive QC reports

## Configuration

Each workflow module can be configured through the main configuration file (`spimquant/config/snakebids.yml`). Key configuration sections include:

- `templates`: Define available templates and atlases
- `registration_level`: Resolution level for registration
- `segmentation_level`: Resolution level for segmentation
- `correction_method`: Intensity correction method (gaussian/n4)
- `seg_method`: Segmentation method (threshold/otsu+k3i2)

## Parallelization

SPIMquant leverages parallelization at multiple levels:

1. **Workflow-level**: Snakemake parallelizes independent rules across subjects
2. **Task-level**: Individual rules use multiple threads for computation
3. **Array-level**: Dask parallelizes array operations within scripts

## Extending the Workflows

To add new functionality:

1. **Add docstrings** to new rules describing their purpose
2. **Follow BIDS naming** conventions for inputs/outputs
3. **Update documentation** when adding significant new features
4. **Test with dry-run** (`-n` flag) before full execution

The modular design makes it easy to add new rules or entire workflow modules while maintaining the existing functionality.
