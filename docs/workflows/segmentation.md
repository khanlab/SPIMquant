# Segmentation and Quantification

The segmentation workflow detects pathology markers (e.g., plaques, cells) and computes quantitative metrics at the object and region level.

## Purpose

Segmentation and quantification enable:
- Automated detection of pathology (amyloid plaques, alpha-synuclein, microglia)
- Quantification of disease burden (counts, density, volume fraction)
- Multi-channel colocalization analysis
- Atlas-based regional statistics
- Group-level comparisons of pathology distribution

The workflow handles multi-channel SPIM data with different stain types and produces comprehensive statistics for downstream analysis.

## Workflow Location

**File**: `spimquant/workflow/rules/segmentation.smk`

## Workflow Stages

### 1. Intensity Correction
Applies bias field correction (Gaussian or N4) to normalize intensities across the field of view. This removes systematic intensity gradients that would affect thresholding.

### 2. Segmentation
Detects objects using threshold-based or multi-Otsu methods. Produces binary masks identifying positive voxels.

### 3. Segmentation Cleaning
Removes edge artifacts and small spurious detections using connected component analysis and spatial filtering.

### 4. Region Properties
Extracts quantitative features for each detected object: size, intensity statistics, centroid location, shape descriptors.

### 5. Coordinate Transformation
Transforms object coordinates from subject to template space using registration warps for cross-subject comparison.

### 6. Colocalization Analysis
Identifies overlapping or nearby objects across channels to detect co-occurring pathologies.

### 7. Atlas-Based Quantification
Maps objects and metrics to atlas regions, computing per-region counts, densities, and volume fractions.

### 8. Visualization
Creates volumetric density maps and subject-level statistical maps for quality control and presentation.

## Key Rule Types

### Correction Rules
- Gaussian bias field correction (fast, approximate)
- N4 bias field correction (slower, more accurate)
- Bias field visualization

### Segmentation Rules
- Simple thresholding
- Multi-Otsu with k-means clustering
- Threshold selection visualization

### Cleaning Rules
- Edge artifact removal
- Connected component filtering
- Size-based object filtering

### Analysis Rules
- Region properties computation
- Spatial coordinate extraction
- Intensity statistics per object

### Transformation Rules
- Point cloud transformation to template space
- Warp field application to coordinates

### Quantification Rules
- Object counts per region
- Density calculation (objects per volume)
- Field fraction (pathology volume fraction)
- Colocalization metrics

### Mapping Rules
- Atlas-based aggregation
- TSV export for statistics
- NIfTI density map generation

## Supported Stains

### Registration Stains
- **PI** (Propidium Iodide): Nuclear stain
- **YOPRO**: Nuclear stain
- **AutoF** (Autofluorescence): Tissue structure

### Segmentation Stains
- **Abeta** / **BetaAmyloid**: Amyloid-beta plaques
- **AlphaSynuclein**: Lewy bodies and neurites
- **Iba1**: Microglia
- **ChAT**: Cholinergic neurons
- **Custom stains**: Configurable in config file

## Segmentation Methods

### Threshold Method

Simple intensity thresholding:

```yaml
seg_method: "threshold"
seg_threshold: 1000  # Intensity value
```

**Use when**:
- Threshold is known from validation
- Background is uniform
- Speed is priority

### Multi-Otsu Method

Automatic threshold selection with k-means:

```yaml
seg_method: "otsu+k3i2"  # k=3 classes, i=2 (use 2nd threshold)
seg_hist_bins: 256
seg_hist_range: [0, 65535]
```

**Use when**:
- Automatic threshold selection needed
- Multiple intensity populations present
- Robust to intensity variations

## Configuration

### Intensity Correction

```yaml
correction_method: "gaussian"  # or "n4"
```

**Gaussian**: Fast, approximate correction using Gaussian blur
**N4**: Accurate correction using ANTs N4 algorithm (slower)

### Segmentation Parameters

```yaml
segmentation_level: 0          # Resolution for segmentation (0=highest)
registration_level: 2          # Resolution for registration
seg_method: "otsu+k3i2"       # Segmentation method
seg_threshold: null            # Manual threshold (null=auto)
```

### Region Property Filters

```yaml
regionprop_filters:
  area:
    min: 10      # Minimum object size (voxels)
    max: 100000  # Maximum object size
  
regionprop_outputs:
  - "label"
  - "area"
  - "centroid"
  - "mean_intensity"
  - "max_intensity"
```

### Colocalization Parameters

```yaml
# In script params (can be exposed to config)
search_radius_multiplier: 1.0  # Scale factor for search distance
overlap_threshold: 0.0         # Minimum overlap to record
```

## Data Types and Formats

### Inputs
- **OME-Zarr** (`.ome.zarr`): Multi-channel SPIM data
- **Corrected Images** (`.ome.zarr`): Bias-corrected intensity
- **Atlas Labels** (`.nii.gz`): Segmentation in subject space
- **Registration Transforms**: Warps for coordinate mapping

### Intermediate Outputs
- **Binary Masks** (`.ome.zarr`): Segmentation results
- **Cleaned Masks** (`.ome.zarr`): After artifact removal
- **Exclusion Masks** (`.ome.zarr`): Removed artifacts

### Final Outputs
- **Region Properties** (`.parquet`): Per-object features
- **Colocalization** (`.parquet`): Cross-channel associations
- **Segmentation Statistics** (`.tsv`): Per-region summaries
- **Density Maps** (`.nii.gz`): Volumetric visualization
- **Field Fraction Maps** (`.nii.gz`): Pathology burden maps

## Metrics and Statistics

### Object-Level Metrics
- **Count**: Number of detected objects
- **Size**: Volume in voxels or µm³
- **Intensity**: Mean, max, integrated intensity
- **Location**: Centroid coordinates (x, y, z)
- **Shape**: Eccentricity, solidity, extent

### Region-Level Metrics
- **Count**: Objects per atlas region
- **Density**: Objects per mm³
- **Field Fraction**: % volume occupied by pathology
- **Total Volume**: Pathology volume in region
- **Mean Intensity**: Average signal intensity

### Colocalization Metrics
- **Overlap**: Spatial overlap between channels
- **Distance**: Nearest neighbor distances
- **Density**: Colocalized objects per region

## Common Workflows

### Standard Segmentation Pipeline

```bash
# Full pipeline with multi-Otsu
spimquant /path/to/data /path/to/output participant \
  --segmentation_level 0 \
  --correction_method gaussian \
  --seg_method otsu+k3i2
```

### Manual Threshold

```bash
# Use predetermined threshold value
spimquant /path/to/data /path/to/output participant \
  --seg_method threshold \
  --seg_threshold 1500
```

### High-Resolution Segmentation

```yaml
# In config
segmentation_level: 0      # Full resolution
registration_level: 2      # Faster registration
```

### Multi-Channel Analysis

```yaml
# Segment multiple stains
stains_for_seg:
  - "Abeta"
  - "Iba1"
  - "AlphaSynuclein"
```

## Quality Control

### Visual Inspection

1. **Threshold plots**: Check histogram and selected threshold
2. **Mask overlays**: Verify segmentation on original image
3. **Cleaned masks**: Ensure artifacts are removed
4. **Density maps**: Check for expected pathology distribution

### Quantitative Checks

- **Object counts**: Should be reasonable for pathology type
- **Size distribution**: Check for spurious small objects
- **Regional variation**: Verify expected anatomical patterns
- **Cross-subject consistency**: Compare metrics across cohort

### Common Issues

**Over-segmentation**:
- Increase threshold
- Adjust Otsu parameters (k, i)
- Improve intensity correction
- Increase minimum object size filter

**Under-segmentation**:
- Decrease threshold
- Check for imaging artifacts
- Verify stain quality
- Adjust intensity correction

**Edge artifacts**:
- Verify cleaning is applied
- Adjust max_extent parameter
- Check registration quality

**Inconsistent results**:
- Standardize intensity correction
- Use consistent resolution levels
- Validate threshold selection
- Check for batch effects

## Performance Considerations

- **Memory**: High-resolution segmentation requires significant RAM (64GB+ recommended)
- **Threads**: Dask parallelization scales well to 128+ threads
- **Storage**: OME-Zarr outputs marked as temp to save disk space
- **Runtime**: ~30-60 minutes per subject depending on resolution and dataset size

## Output File Organization

```
derivatives/
├── sub-001/
│   ├── micr/
│   │   ├── sub-001_stain-Abeta_desc-otsu+k3i2+cleaned_regionprops.parquet
│   │   ├── sub-001_stain-Abeta_level-2_desc-otsu+k3i2+cleaned_fieldfrac.nii.gz
│   │   ├── sub-001_seg-AllenCCFv3_from-ABAv3_desc-otsu+k3i2+cleaned_mergedsegstats.tsv
│   │   └── sub-001_space-ABAv3_desc-otsu+k3i2+cleaned_coloc.parquet
```

## Dependencies

- **Dask**: Parallel array processing
- **scikit-image**: Segmentation algorithms (Otsu, connected components)
- **pandas**: Data table handling
- **nibabel**: NIfTI I/O
- **zarrnii**: OME-Zarr handling
- **ANTs**: N4 correction (optional)

## Related Workflows

- **[Import](import.md)**: Provides input data
- **[Masking](masking.md)**: Creates tissue masks
- **[Registration](templatereg.md)**: Provides transforms and atlas labels
- **[Group Statistics](groupstats.md)**: Analyzes segmentation statistics across subjects
- **[Patches](patches.md)**: Extracts patches from segmentation results
