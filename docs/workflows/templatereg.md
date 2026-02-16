# Template Registration

The template registration workflow aligns subject SPIM data to a standard reference template using affine and deformable registration.

## Purpose

Template registration enables:
- Spatial normalization across subjects for group analysis
- Atlas-based anatomical labeling
- Comparison of pathology distributions across cohorts
- Population-level statistical mapping
- Cross-study data integration

The workflow implements a multi-stage registration pipeline optimized for SPIM microscopy data, handling the large data volumes and intensity characteristics specific to cleared tissue imaging.

## Workflow Location

**File**: `spimquant/workflow/rules/templatereg.smk`

## Workflow Stages

### 1. Intensity Correction
Applies N4 bias field correction to remove intensity inhomogeneities that would degrade registration accuracy. The bias field is estimated within the brain mask and applied to create a more uniform intensity distribution.

### 2. Affine Registration
Performs global alignment using 12-parameter affine transformation (translation, rotation, scaling, shearing). Uses normalized mutual information (NMI) as the similarity metric, which is robust to intensity differences between modalities.

### 3. Deformable Registration
Refines alignment with non-linear deformation field to capture local anatomical variations. Uses multi-resolution optimization with Gaussian smoothing to avoid local minima and achieve stable convergence.

### 4. Label Transformation
Applies registration transforms to move atlas segmentations between template and subject space. Forward transforms map subject data to template; inverse transforms bring atlas labels to subject space.

### 5. Quality Control
Generates interactive HTML report with checkerboard visualizations, edge overlays, and warp field visualizations to assess registration quality.

## Key Rule Types

### Preprocessing Rules
- N4 bias field correction
- Brain mask application
- Template cropping (optional, for hemisphere-specific analysis)

### Registration Rules
- Affine registration (greedy, 12-DOF)
- Deformable registration (greedy, multi-resolution)
- Transform format conversion (RAS ↔ ITK)

### Transformation Rules
- Apply transforms to images (linear interpolation)
- Apply transforms to labels (nearest-neighbor interpolation)
- Warp OME-Zarr to template space
- Resample at multiple resolutions

### Quality Control Rules
- Registration QC report generation
- Visual overlay creation

## Supported Templates

SPIMquant supports multiple standard templates:

- **ABAv3**: Allen Brain Atlas version 3 (25 μm mouse)
- **gubra**: Gubra LSFM template (20 μm mouse)
- **MBMv3**: Marmoset Brain Mapping v3 (150 μm marmoset)
- **turone**: Ontario Human Brain Atlas
- **MouseIn**: MRI mouse template

Each template includes:
- Anatomical reference image
- Brain mask
- One or more atlas parcellations
- Label lookup tables

## Configuration

### Template Selection

```yaml
template: "ABAv3"  # Primary template for registration
```

### Registration Parameters

```yaml
registration_level: 2        # Resolution level for registration
templatereg:
  desc: "N4brain"           # Use N4-corrected, masked image
  stain: "PI"               # Stain for registration (e.g., PI, YOPRO, AutoF)
```

### Hemisphere Cropping

```yaml
template_crop: "left"       # Register to left hemisphere only
# or "right" for right hemisphere
# or null for whole brain
```

### Performance Parameters

```yaml
sloppy: false               # Set true for quick testing (reduced iterations)
```

When `sloppy: true`:
- Affine iterations: 10x0x0 (vs 100x100)
- Deformable iterations: 10x0x0 (vs 100x50)

## Data Types and Formats

### Inputs
- **Subject SPIM** (`.nii.gz`, `.ome.zarr`): Multi-channel microscopy data
- **Template Anatomy** (`.nii.gz`): Reference image
- **Template Mask** (`.nii.gz`): Template brain mask
- **Atlas Segmentations** (`.nii.gz`): Parcellation files

### Outputs

**Transformation Files**:
- Affine transform (`.txt`): RAS and ITK formats
- Deformation field (`.nii.gz`): Forward and inverse warps

**Warped Images**:
- Affine-warped subject (`.nii.gz`)
- Deformable-warped subject (`.nii.gz`, `.ome.zarr`)
- Atlas in subject space (`.nii.gz`)

**Quality Control**:
- Registration QC report (`.html`)

## Algorithm Details

### Affine Registration

Uses greedy's affine registration:
- **Initialization**: Image center alignment
- **Metric**: NMI (handles intensity differences)
- **Optimization**: Multi-resolution (100x100 iterations)
- **DOF**: 12 (full affine)

### Deformable Registration

Uses greedy's deformable registration:
- **Metric**: NMI (configurable)
- **Smoothing**: Multi-scale (σ=4vox, σ=2vox)
- **Iterations**: 100x50 (coarse-to-fine)
- **Regularization**: Gaussian smoothing of deformation field

### Transform Composition

Transforms are composed for multi-stage warping:

```
Template → Subject: [Inverse Warp] → [Inverse Affine]
Subject → Template: [Affine] → [Forward Warp]
```

## Multi-Resolution Strategy

Registration uses coarser resolution for efficiency:

1. **Downsampled registration** (level 2-3):
   - Faster computation
   - More robust to noise
   - Avoids local minima

2. **Full-resolution application**:
   - Transforms are resolution-independent
   - Can warp level 0 data using transforms from level 2
   - Balances accuracy and speed

## OME-Zarr Handling

The workflow handles OME-Zarr formats efficiently:

- **Lazy loading**: Only loads necessary resolution levels
- **Chunked processing**: Uses Dask for out-of-core computation
- **Multi-resolution output**: Generates OME-Zarr pyramids in template space

## Common Workflows

### Standard Registration

```bash
# Full-resolution registration
spimquant /path/to/data /path/to/output participant \
  --template ABAv3 \
  --registration_level 2
```

### Fast Test Run

```bash
# Quick test with reduced iterations
spimquant /path/to/data /path/to/output participant \
  --template ABAv3 \
  --sloppy
```

### Hemisphere-Specific Analysis

```bash
# Register left hemisphere only
spimquant /path/to/data /path/to/output participant \
  --template ABAv3 \
  --template_crop left
```

### Custom Stain Selection

```yaml
# In configuration file
templatereg:
  stain: "YOPRO"  # Use YOPRO channel for registration
```

## Quality Control

### Visual Assessment

The registration QC report includes:
1. **Checkerboard comparison**: Interleaved template/subject patches
2. **Edge overlay**: Template edges on subject image
3. **Warp field visualization**: Deformation magnitude and direction
4. **Atlas overlay**: Template labels on warped subject

### Quantitative Metrics

Check for:
- Proper alignment of major structures (cortex, ventricles, white matter)
- Smooth deformation fields without folding
- Consistent correspondence across subjects

### Common Issues

**Poor affine alignment**:
- Check intensity normalization
- Verify correct stain selection
- Ensure brain mask quality
- Try different initialization

**Deformable registration artifacts**:
- Reduce regularization (increase smoothing)
- Adjust iterations
- Check for intensity outliers
- Verify affine quality first

**Inconsistent results**:
- Standardize preprocessing (N4, masking)
- Use consistent resolution levels
- Verify template appropriateness

## Performance Considerations

- **Memory**: Deformable registration is memory-intensive (~16GB+)
- **Threads**: Linear speedup up to 32 threads
- **Runtime**: ~5-15 minutes per subject (depending on resolution)
- **Storage**: OME-Zarr outputs can be large; use temp directories

## Dependencies

- **greedy**: Affine and deformable registration
- **ANTs**: N4 correction, transform application
- **c3d**: Transform format conversion
- **zarrnii**: OME-Zarr handling
- **Dask**: Parallel array processing

## Related Workflows

- **[Import](import.md)**: Provides templates and atlases
- **[Masking](masking.md)**: Creates brain masks for registration
- **[Segmentation](segmentation.md)**: Uses registration for atlas-based analysis
- **[Group Statistics](groupstats.md)**: Analyzes data in common template space
