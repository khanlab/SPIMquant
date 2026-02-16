# Brain Masking

The masking workflow creates brain tissue masks by combining Gaussian mixture modeling with template-derived spatial priors.

## Purpose

Brain masking is critical for:
- Removing background and non-brain tissue
- Improving registration accuracy
- Restricting analysis to relevant tissue regions
- Reducing computational requirements

The workflow uses a hybrid approach combining intensity-based tissue classification (Atropos) with template-based spatial priors for robust brain extraction.

## Workflow Location

**File**: `spimquant/workflow/rules/masking.smk`

## Workflow Stages

### 1. Image Preprocessing
Downsamples and normalizes the SPIM image for efficient tissue classification. Applies log transformation to improve separation of intensity classes.

### 2. Tissue Classification
Uses ANTs Atropos to classify voxels into k intensity classes using Gaussian mixture modeling (GMM) with Markov random field (MRF) spatial smoothing. The MRF encourages spatial coherence in the classification.

### 3. Initial Registration
Performs fast affine registration to the template to obtain a rough spatial alignment. This enables template masks to be used as anatomical priors.

### 4. Template Prior Transformation
Warps the template brain mask to subject space using the inverse affine transformation. This provides anatomical prior knowledge about expected brain location and shape.

### 5. Mask Combination
Combines GMM-based tissue classes with the template prior to create the final brain mask. This hybrid approach leverages both intensity information and anatomical priors.

## Key Rule Types

### Preprocessing Rules
- Downsampling for computational efficiency
- Intensity normalization and transformation
- Initial foreground mask creation

### Classification Rules
- K-class Gaussian mixture modeling
- MRF spatial smoothing
- Posterior probability map generation

### Registration Rules
- Affine template-to-subject alignment
- Transform template masks to subject space
- Inverse transformation handling

### Mask Creation Rules
- GMM class combination
- Prior-guided refinement
- Alternative masking strategies (GMM-only)

## Configuration Parameters

### GMM Parameters

```yaml
masking:
  gmm_k: 3                    # Number of intensity classes
  gmm_bg_class: 1             # Background class label
  pre_atropos_downsampling: "10%"  # Downsampling factor
```

### MRF Parameters

MRF smoothing parameters are defined in the rule:
- `mrf_smoothing`: Controls strength of spatial regularization (default: 0.3)
- `mrf_radius`: Neighborhood size for MRF (default: 2x2x2 voxels)

## Registration Parameters

Initial affine registration uses:
- **DOF**: 12 (full affine)
- **Metric**: NMI (Normalized Mutual Information)
- **Iterations**: Configurable via `sloppy` mode

## Data Types and Formats

### Inputs
- **SPIM NIfTI** (`.nii.gz`): Subject SPIM image at registration level
- **Template** (`.nii.gz`): Reference template anatomy
- **Template Mask** (`.nii.gz`): Template brain mask

### Intermediate Outputs
- **Downsampled Image** (`.nii.gz`): Preprocessed for Atropos
- **Segmentation** (`.nii`): K-class tissue classification
- **Posteriors** (`.nii`): Probability maps per class
- **Affine Transform** (`.txt`): Subject-to-template alignment
- **Template Mask in Subject Space** (`.nii.gz`): Warped prior

### Final Outputs
- **Brain Mask** (`.nii.gz`): Binary brain tissue mask

## Algorithm Details

### Gaussian Mixture Modeling

Atropos fits a GMM to the image intensities:

```
P(I|class_k) = N(μ_k, σ_k²)
```

Where:
- `I`: Voxel intensity
- `μ_k`: Mean intensity of class k
- `σ_k`: Standard deviation of class k

### MRF Spatial Prior

The MRF adds spatial regularization:

```
P(class_k|neighbors) ∝ exp(-β * Σ(class_i ≠ class_k))
```

This encourages neighboring voxels to have the same class label.

### Prior Integration

The final mask combines GMM classes with the template prior:

```python
brain_mask = (gmm_classes ∈ brain_classes) ∧ (template_mask > threshold)
```

## Common Workflows

### Standard Masking

```bash
# Default 3-class GMM with template prior
spimquant /path/to/data /path/to/output participant
```

### Fast Masking (Sloppy Mode)

```bash
# Reduced iterations for quick testing
spimquant /path/to/data /path/to/output participant --sloppy
```

### Tuning GMM Parameters

Adjust parameters in configuration:

```yaml
masking:
  gmm_k: 4              # Try 4 classes for complex intensity distributions
  gmm_bg_class: 1       # Specify which class is background
```

## Quality Control

### Visual Inspection

Check mask quality by overlaying on original image:
1. Load subject SPIM image
2. Load brain mask as overlay
3. Verify brain boundaries are captured
4. Check for inclusion of non-brain tissue

### Common Issues

**Mask too restrictive**:
- Increase number of GMM classes
- Adjust background class label
- Check template mask accuracy

**Mask includes non-brain tissue**:
- Improve initial registration
- Use more restrictive template mask
- Adjust MRF smoothing strength

**Inconsistent results across subjects**:
- Verify consistent stain selection
- Check for intensity normalization issues
- Ensure template is appropriate for the data

## Dependencies

- **ANTs**: Atropos tissue classification, N4 correction
- **greedy**: Affine registration
- **c3d**: Image manipulation and downsampling
- **nibabel**: NIfTI I/O

## Performance Considerations

- **Downsampling**: 10% downsampling provides good balance of speed and accuracy
- **MRF radius**: Larger radius increases smoothing but also computation time
- **Number of classes**: More classes = longer computation but better separation

## Related Workflows

- **[Import](import.md)**: Provides template masks
- **[Registration](templatereg.md)**: Uses masks for intensity correction
- **[Segmentation](segmentation.md)**: May use masks to restrict segmentation
