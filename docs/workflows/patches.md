# Patch Extraction

The patch extraction workflow creates fixed-size 3D volumes from SPIM data and segmentation masks based on atlas region sampling.

## Purpose

Patch extraction enables:
- Training machine learning models on region-specific data
- Visual quality control of segmentation results
- Creating datasets for downstream analysis
- Extracting high-resolution regions of interest
- Generating examples for presentations and publications

Patches are sampled from specific atlas regions with reproducible randomization, making them ideal for systematic analysis and validation.

## Workflow Location

**File**: `spimquant/workflow/rules/patches.smk`

**Note**: This file contains comprehensive docstrings for all rules describing their functionality.

## Workflow Types

The workflow provides four types of patch extraction:

### 1. Raw SPIM Patches
Extracts patches from original, uncorrected SPIM data at the highest resolution (level 0).

**Use for**:
- Visualizing raw data quality
- Comparing preprocessing methods
- Training models on raw intensities

### 2. Corrected SPIM Patches
Extracts patches from intensity-corrected SPIM data (Gaussian or N4).

**Use for**:
- Assessing correction quality
- Training on normalized data
- Comparing correction methods

### 3. Segmentation Mask Patches
Extracts binary mask patches from cleaned segmentation results.

**Use for**:
- Validating segmentation accuracy
- Training segmentation models
- Creating ground truth datasets

### 4. Imaris Crops
Creates high-resolution Imaris-format datasets for regions of interest.

**Use for**:
- Interactive 3D visualization in Imaris
- Detailed examination of specific regions
- Creating publication-quality renderings

## Configuration

### Patch Parameters

```yaml
patch_size: [256, 256, 256]       # Size of patches in voxels [x, y, z]
n_patches_per_label: 10           # Number of patches to sample per region
patch_seed: 42                     # Random seed for reproducibility
```

### Region Selection

```yaml
patch_labels:                      # Atlas regions to sample from
  - "Hippocampus"
  - "Cortex"
  - "Striatum"
  # or null to sample from all regions
```

### Imaris Crop Parameters

```yaml
crop_labels:                       # Regions for Imaris export
  - "CA1"
  - "DG"
  # Typically fewer, larger regions for detailed visualization
```

## Data Types and Formats

### Inputs
- **OME-Zarr** (`.ome.zarr`): Multi-resolution SPIM data
- **Segmentation Masks** (`.ome.zarr`): Binary detection masks
- **Atlas Labels** (`.nii.gz`): Region definitions in subject space
- **Label Table** (`.tsv`): Region names and abbreviations

### Outputs
- **Patch Directories**: Collections of NIfTI patches
  - Named: `{atlas}_{label_abbrev}_patch{n}.nii.gz`
- **Imaris Datasets**: Hierarchical multi-resolution formats
  - Organized by region in directories

## Sampling Strategy

### Random Sampling Within Regions

For each specified atlas region:
1. Identify all voxels belonging to that region
2. Randomly sample N center points
3. Extract fixed-size patches around each center
4. Ensure patches don't extend beyond image boundaries
5. Save with descriptive filenames

### Reproducibility

The `patch_seed` parameter ensures:
- Consistent sampling across runs
- Reproducible train/test splits
- Comparable validation results

### Spatial Distribution

Sampling is uniform within each region but respects:
- Region boundaries (patches don't cross regions)
- Image boundaries (no edge artifacts)
- Resolution level (appropriate for analysis type)

## Resolution Levels

### High-Resolution Patches (Level 0)
- Maximum detail for object detection
- Large file sizes
- Best for training segmentation models

### Downsampled Patches (Level 1-2)
- Faster to process
- Smaller file sizes
- Sufficient for intensity-based analysis

## Output Organization

```
derivatives/
└── sub-001/
    └── micr/
        ├── sub-001_stain-Abeta_seg-AllenCCFv3_from-ABAv3_level-0_desc-raw_SPIM.patches/
        │   ├── AllenCCFv3_CA1_patch001.nii.gz
        │   ├── AllenCCFv3_CA1_patch002.nii.gz
        │   ├── AllenCCFv3_CTX_patch001.nii.gz
        │   └── ...
        ├── sub-001_stain-Abeta_seg-AllenCCFv3_from-ABAv3_level-0_desc-otsu+k3i2+cleaned_mask.patches/
        │   ├── AllenCCFv3_CA1_patch001.nii.gz
        │   └── ...
        └── sub-001_seg-AllenCCFv3_from-ABAv3_level-0_desc-crop_SPIM.imaris/
            ├── CA1/
            └── DG/
```

## Common Workflows

### Extract Patches for Validation

```yaml
# Configuration
patch_labels:
  - "Hippocampus"
  - "Cortex"
n_patches_per_label: 20
patch_size: [256, 256, 256]
```

```bash
# Run workflow
spimquant /path/to/data /path/to/output participant \
  --segmentation_level 0
```

### Create Imaris Datasets

```yaml
# Configuration
crop_labels:
  - "CA1"          # Hippocampal CA1
  - "DG"           # Dentate gyrus
```

Imaris crops include full resolution and create optimized hierarchical storage for interactive visualization.

### Matched Raw and Mask Patches

Extract both raw and segmentation patches with matching sampling:

```yaml
patch_seed: 42    # Same seed for both
n_patches_per_label: 50
```

This creates paired datasets for model training.

## Quality Control

### Visual Inspection

1. Load patches in image viewer (e.g., ITK-SNAP, FIJI)
2. Verify correct region sampling
3. Check patch quality (contrast, artifacts)
4. Validate mask accuracy

### Systematic Validation

```python
# Example Python code to validate patches
import nibabel as nib
from pathlib import Path

patch_dir = Path("derivatives/sub-001/.../patches/")
for patch_file in patch_dir.glob("*.nii.gz"):
    img = nib.load(patch_file)
    data = img.get_fdata()
    
    # Check patch size
    assert data.shape == (256, 256, 256)
    
    # Check intensity range
    print(f"{patch_file.name}: min={data.min():.1f}, max={data.max():.1f}")
```

### Common Issues

**Insufficient sampling locations**:
- Region too small for requested number of patches
- Solution: Reduce n_patches_per_label or increase region size

**Patches at image edges**:
- Sampling too close to boundaries
- Solution: Script automatically handles this, but verify visually

**Mismatched patches between modalities**:
- Different random seeds
- Solution: Ensure consistent patch_seed across runs

## Use Cases

### Machine Learning

**Segmentation Model Training**:
1. Extract raw SPIM patches
2. Extract corresponding mask patches
3. Use as input/target pairs for U-Net, etc.

**Classification**:
1. Extract region-specific patches
2. Label by region or pathology type
3. Train classifier (CNN, ResNet)

### Quality Control

**Segmentation Validation**:
1. Extract patches with high/low confidence detections
2. Manual review by expert
3. Adjust parameters based on findings

**Cross-Subject Comparison**:
1. Extract patches from same region across subjects
2. Visual comparison of data quality
3. Identify batch effects or outliers

### Visualization

**Publication Figures**:
1. Extract representative patches
2. Create montages or galleries
3. Overlay masks on raw data

**Presentations**:
1. Export Imaris crops
2. Create 3D renderings
3. Generate rotation videos

## Performance Considerations

- **Memory**: Patches are extracted sequentially (low memory)
- **Storage**: 256³ patches at level 0 can be large (~134 MB each)
- **Runtime**: ~30-60 minutes depending on number of patches
- **Parallelization**: Extraction uses multi-threading efficiently

## Dependencies

- **zarrnii**: OME-Zarr reading and patch extraction
- **nibabel**: NIfTI I/O
- **numpy**: Array operations
- **dask**: Lazy loading for memory efficiency
- **scikit-image**: Label region processing

## Related Workflows

- **[Segmentation](segmentation.md)**: Produces masks for patch extraction
- **[Registration](templatereg.md)**: Provides atlas labels for region sampling
- **[Import](import.md)**: Provides label definitions

## Advanced Topics

### Custom Sampling Strategies

Modify the sampling script to implement:
- Density-weighted sampling (more patches in high-density regions)
- Stratified sampling (balanced across subjects)
- Active learning sampling (patches with high uncertainty)

### Integration with External Tools

Patches can be directly loaded in:
- **PyTorch**: `nibabel` → `numpy` → `torch.Tensor`
- **TensorFlow**: `nibabel` → `numpy` → `tf.Tensor`
- **napari**: Interactive visualization and annotation
- **FIJI/ImageJ**: Batch processing and analysis

### Patch Datasets for Sharing

Package extracted patches for data sharing:
```bash
# Create dataset archive
tar -czf patches_dataset.tar.gz derivatives/*/micr/*.patches/

# Include metadata
cp participants.tsv patches_dataset/
cp label_table.tsv patches_dataset/
```
