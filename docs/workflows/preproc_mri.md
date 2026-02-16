# MRI Preprocessing and Cross-Modality Registration

The MRI preprocessing workflow handles co-registration of in-vivo MRI with ex-vivo SPIM data for multi-modal analysis.

## Purpose

MRI-SPIM co-registration enables:
- Assessment of tissue changes due to perfusion fixation
- Quantification of tissue shrinkage or expansion
- Multi-modal anatomical validation
- Complementary structural information
- Longitudinal in-vivo to ex-vivo tracking

This workflow is optional and used only when both MRI and SPIM data are available for the same subject.

## Workflow Location

**File**: `spimquant/workflow/rules/preproc_mri.smk`

**Note**: This file contains docstrings for key rules describing their functionality.

## Workflow Stages

### 1. MRI Preprocessing
Applies N4 bias field correction to T2-weighted MRI to normalize intensities and improve subsequent registration.

### 2. MRI to Template Registration
Registers MRI to a standard MRI template (e.g., MouseIn) to obtain brain masks via template propagation.

### 3. Brain Mask Transformation
Warps template brain mask to MRI space for brain extraction.

### 4. MRI Brain Extraction
Applies brain mask to create skull-stripped MRI.

### 5. MRI to SPIM Registration
Registers brain-extracted MRI to SPIM using rigid or affine followed by deformable registration.

### 6. Parameter Optimization
Provides rules for testing multiple parameter combinations to optimize registration quality.

### 7. Transformation Concatenation
Chains transformations to enable direct MRI-to-template warping via SPIM space.

### 8. Jacobian Analysis
Computes Jacobian determinant of the deformation field to quantify local tissue volume changes.

## Key Rule Types

### Preprocessing Rules
- N4 bias field correction for MRI
- Intensity normalization
- Brain mask creation

### Registration Rules
- MRI to MRI template (rigid + deformable)
- MRI to SPIM (rigid/affine + deformable)
- Template mask to MRI transformation

### Parameter Tuning Rules
- Grid search over registration parameters
- Multiple iteration schemes
- Various smoothing kernels and radii

### Transformation Rules
- Concatenated transformations (MRI → SPIM → Template)
- Brain mask warping with Jacobian computation
- Multi-stage transform composition

### Quality Control Rules
- Visual comparison of registration results
- Jacobian determinant maps
- Overlap metrics

## Supported MRI Templates

### MouseIn Template
Standard MRI template for mouse brain:
- **Resolution**: ~100 μm isotropic
- **Modality**: T2-weighted
- **Format**: NIfTI
- **Includes**: Brain mask, tissue segmentation

## Registration Strategy

### Multi-Stage Approach

1. **MRI Template Registration**:
   - Purpose: Obtain brain mask prior
   - Method: Rigid (6-DOF) + deformable
   - Metric: NCC (Normalized Cross-Correlation)

2. **MRI to SPIM Registration**:
   - Purpose: Align modalities
   - Method: Affine (12-DOF) + deformable
   - Metric: NMI (handles intensity differences) + NCC

3. **Concatenated Warping**:
   - Purpose: Direct MRI to atlas transformation
   - Method: Compose all transforms
   - Benefit: Single resampling step

## Configuration

### Template Selection

```yaml
template_mri: "MouseIn"     # MRI template for brain masking
template: "ABAv3"           # SPIM template for final alignment
```

### Registration Parameters

The workflow uses wildcards for parameter exploration:

```yaml
# These are specified as wildcards in filenames
iters: "100x100x50x0"       # Multi-resolution iterations
dof: "12"                   # Degrees of freedom (6=rigid, 12=affine)
radius: "2x2x2"             # NCC radius
gradsigma: "3"              # Gradient smoothing
warpsigma: "3"              # Warp smoothing
```

### Parameter Tuning

Enable comprehensive parameter search:

```yaml
# Create all combinations
iters: ["100x50x50", "100x100x50x0"]
dof: [6, 12]
radius: ["2x2x2", "3x3x3", "4x4x4"]
gradsigma: [2, 3, 4, 5]
warpsigma: [2, 3, 4, 5]
```

## Data Types and Formats

### Inputs
- **T2w MRI** (`.nii.gz`): In-vivo structural MRI
- **SPIM Data** (`.nii.gz`): Ex-vivo lightsheet data
- **MRI Template** (`.nii.gz`): Reference MRI
- **Template Mask** (`.nii.gz`): Brain mask

### Intermediate Outputs
- **N4-Corrected MRI** (`.nii.gz`)
- **MRI Brain Mask** (`.nii.gz`)
- **Affine Transforms** (`.txt`)
- **Deformation Fields** (`.nii.gz`)

### Final Outputs
- **Brain-Extracted MRI** (`.nii.gz`)
- **MRI in SPIM Space** (`.nii.gz`)
- **MRI in Template Space** (`.nii.gz`)
- **Jacobian Map** (`.nii.gz`): Volume change map

## Jacobian Determinant Interpretation

The Jacobian determinant quantifies local volume changes:

- **J > 1**: Tissue expansion (e.g., due to perfusion, swelling)
- **J = 1**: No volume change (perfect preservation)
- **J < 1**: Tissue shrinkage (e.g., due to dehydration, clearing)

Common observations:
- **Cleared tissue**: Often shows J < 1 (10-30% shrinkage)
- **Gray matter**: More shrinkage than white matter
- **Ventricles**: May show expansion if filled during perfusion

## Common Workflows

### Basic MRI-SPIM Co-Registration

```yaml
# Configuration
template_mri: "MouseIn"
template: "ABAv3"
```

```bash
# Process subject with both MRI and SPIM
spimquant /path/to/data /path/to/output participant \
  --filter-T2w acquisition=highres
```

### Parameter Optimization

```bash
# Generate all parameter combinations for tuning
spimquant /path/to/data /path/to/output participant \
  --targets_by_analysis_level participant=all_tune_mri_spim_reg
```

This creates dozens of registration results for visual comparison.

### Jacobian Analysis

After successful registration, examine tissue deformation:

```python
import nibabel as nib
import numpy as np

# Load Jacobian map
jac_img = nib.load("derivatives/.../desc-brain_jacobian.nii.gz")
jac_data = jac_img.get_fdata()

# Compute statistics
mean_jac = np.mean(jac_data[jac_data > 0])
print(f"Mean volume change: {(mean_jac - 1) * 100:.1f}%")

# Shrinkage: mean_jac < 1
# Expansion: mean_jac > 1
```

## Quality Control

### Visual Assessment

1. **MRI to template alignment**: Check cortex, ventricles, subcortical structures
2. **MRI to SPIM alignment**: Verify correspondence of major features
3. **Jacobian map**: Inspect for smooth gradients (no sharp discontinuities)

### Quantitative Metrics

- **Overlap**: Dice coefficient of brain masks
- **Intensity correlation**: After registration
- **Jacobian statistics**: Mean, std, range
- **Deformation smoothness**: No folding or tears

### Common Issues

**Poor MRI-template alignment**:
- MRI quality issues (motion, artifacts)
- Wrong template selected
- Insufficient intensity contrast
- Try different parameter combinations

**MRI-SPIM misalignment**:
- Large deformation due to tissue processing
- Different anatomical coverage
- Intensity differences too large
- Adjust deformable registration parameters

**Unrealistic Jacobian values**:
- Registration failed (check alignment first)
- Inappropriate smoothing parameters
- Try more regularization (increase warpsigma)

**Parameter tuning takes too long**:
- Reduce parameter grid
- Use sloppy mode for initial screening
- Focus on most important parameters (iters, radius)

## Performance Considerations

- **Memory**: MRI registration less demanding than SPIM (~8-16 GB)
- **Runtime**: 
  - MRI template registration: ~5-15 minutes
  - MRI to SPIM registration: ~10-20 minutes
  - Parameter tuning: hours to days
- **Storage**: MRI files relatively small (<1 GB per subject)

## Biological Applications

### Tissue Shrinkage Studies

Quantify effects of clearing protocols:
```python
# Compare Jacobian between protocols
protocol_a_shrinkage = 1 - mean_jacobian_a
protocol_b_shrinkage = 1 - mean_jacobian_b

print(f"Protocol A: {protocol_a_shrinkage*100:.1f}% shrinkage")
print(f"Protocol B: {protocol_b_shrinkage*100:.1f}% shrinkage")
```

### Longitudinal Studies

Track disease progression with in-vivo MRI followed by ex-vivo SPIM at endpoint:
1. Serial MRI during life
2. Terminal SPIM with pathology staining
3. Co-register all timepoints
4. Correlate MRI changes with histology

### Multi-Modal Validation

Validate SPIM findings with MRI:
- Lesion location correspondence
- Volume measurements
- Structural abnormalities

## Dependencies

- **ANTs**: N4 correction
- **greedy**: Registration (rigid, affine, deformable)
- **c3d**: Transform manipulation, Jacobian calculation
- **nibabel**: NIfTI I/O

## Related Workflows

- **[Import](import.md)**: Provides MRI templates
- **[Masking](masking.md)**: Similar brain masking strategy
- **[Registration](templatereg.md)**: SPIM to template registration
- **[Segmentation](segmentation.md)**: Can leverage MRI-derived anatomical priors

## Advanced Topics

### Custom MRI Templates

Add your own MRI template:

```yaml
templates:
  CustomMRI:
    anat: "path/to/custom_T2w.nii.gz"
    mask: "path/to/custom_mask.nii.gz"
```

### Multi-Contrast MRI

If multiple MRI contrasts available (T1, T2, etc.):
1. Register all to same space
2. Use contrast with best tissue differentiation for SPIM registration
3. Warp all contrasts to SPIM space for multi-modal analysis

### Symmetric Registration

For unbiased volume analysis:
1. Register MRI → SPIM
2. Register SPIM → MRI
3. Compose inverse transforms
4. Compare Jacobians in both directions
