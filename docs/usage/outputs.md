# Output and Results

SPIMquant generates comprehensive outputs in BIDS-derivatives format, including registered images, quantitative measurements, and quality control data. This section describes the output structure and how to interpret results.

## Output Directory Structure

### BIDS-Derivatives Organization

```
output/
├── spimquant/                          # Main derivatives folder
│   ├── dataset_description.json        # Dataset metadata
│   ├── sub-001/                        # Subject-specific results
│   │   └── micr/                       # Microscopy data type
│   │       ├── *_space-{template}_SPIM.nii           # Registered images
│   │       ├── *_seg-{atlas}_fieldfrac.nii          # Field fraction maps
│   │       ├── *_stain-{stain}_blobs.npy            # Detected features
│   │       ├── *_seg-{atlas}_blobdensity.nii        # Density maps
│   │       └── *_desc-{method}_mask.nii             # Brain masks
│   ├── tpl-{template}/                 # Template space files
│   │   ├── tpl-{template}_anat.nii.gz               # Template anatomy
│   │   ├── tpl-{template}_seg-{atlas}_dseg.nii.gz   # Atlas parcellation
│   │   └── tpl-{template}_seg-{atlas}_dseg.tsv      # Region labels
│   └── logs/                           # Processing logs
└── work/                               # Temporary files (can be deleted)
```

## Core Output Files

### Registered Images

**Filename Pattern:**
```
sub-{subject}_sample-{sample}_stain-{stain}_space-{template}_SPIM.nii
```

**Example:**
```
sub-001_sample-brain_stain-PI_space-ABAv3_SPIM.nii
sub-001_sample-brain_stain-abeta_space-ABAv3_SPIM.nii
```

**Description:**
- Original SPIM data registered to template space
- Intensity-corrected and normalized
- Same voxel size as template (typically 25μm isotropic)
- Multiple stains per subject registered to same space

**Usage:**
```bash
# Visualize with ITK-SNAP
itksnap output/spimquant/sub-001/micr/sub-001_sample-brain_stain-PI_space-ABAv3_SPIM.nii

# Overlay with template
itksnap \
  -g output/spimquant/tpl-ABAv3/tpl-ABAv3_anat.nii.gz \
  -o output/spimquant/sub-001/micr/sub-001_sample-brain_stain-abeta_space-ABAv3_SPIM.nii
```

### Field Fraction Maps

**Filename Pattern:**
```
sub-{subject}_seg-{atlas}_stain-{stain}_fieldfrac.nii
sub-{subject}_space-{template}_seg-{atlas}_stain-{stain}_fieldfrac.nii
```

**Example:**
```
sub-001_seg-roi82_stain-abeta_fieldfrac.nii
sub-001_space-ABAv3_seg-roi82_stain-abeta_fieldfrac.nii
```

**Description:**
- Voxel-wise proportion of positive signal (0.0 to 1.0)
- Generated using specified segmentation method
- Available at multiple resolution levels
- Per-region and per-voxel quantification

**Interpretation:**
```
Value 0.0:    No positive signal detected
Value 0.5:    50% of voxel contains positive signal  
Value 1.0:    Entire voxel is positive signal
```

### Blob Detection Results

**Filename Pattern:**
```
sub-{subject}_stain-{stain}_blobs.npy
```

**Example:**
```
sub-001_stain-abeta_blobs.npy
sub-001_stain-Iba1_blobs.npy
```

**Description:**
- Individual detected features (blobs/objects)
- NumPy array format with feature properties
- Coordinates in template space
- Size, intensity, and shape characteristics

**Data Structure:**
```python
import numpy as np
blobs = np.load('sub-001_stain-abeta_blobs.npy')

# Array columns:
# [y, x, z, sigma_y, sigma_x, sigma_z, intensity]
print(f"Number of blobs: {len(blobs)}")
print(f"Average intensity: {np.mean(blobs[:, -1])}")
```

### Blob Density Maps

**Filename Pattern:**
```
sub-{subject}_seg-{atlas}_stain-{stain}_blobdensity.nii
```

**Example:**
```
sub-001_seg-roi82_stain-abeta_blobdensity.nii
```

**Description:**
- Spatial density of detected features
- Units: blobs per cubic millimeter
- Smoothed density estimates
- Atlas-region specific measurements

### Brain Masks

**Filename Pattern:**
```
sub-{subject}_desc-{method}_mask.nii
```

**Example:**
```
sub-001_desc-atropos_mask.nii
sub-001_desc-brain_mask.nii
```

**Description:**
- Binary masks defining brain tissue
- Generated using tissue segmentation
- Used for constraining analysis regions
- Quality control for registration

## Template Space Files

### Template Anatomy

**Filename:**
```
tpl-{template}_anat.nii.gz
```

**Description:**
- Reference anatomy for visualization
- T2-weighted or equivalent contrast
- Standard space for all registered data
- Used as background for overlays

### Atlas Parcellations

**Filename:**
```
tpl-{template}_seg-{atlas}_dseg.nii.gz
```

**Description:**
- Integer-valued brain region labels
- Different resolutions available (roi22, roi82, roi198)
- Corresponds to label table (.tsv file)
- Used for region-of-interest analysis

### Label Tables

**Filename:**
```
tpl-{template}_seg-{atlas}_dseg.tsv
```

**Format:**
```tsv
index	name	abbreviation	color
1	Cortex	Ctx	255,0,0
2	Hippocampus, CA1	CA1	0,255,0
3	Hippocampus, CA3	CA3	0,128,255
```

**Usage:**
```python
import pandas as pd

# Load label table
labels = pd.read_csv('tpl-ABAv3_seg-roi82_dseg.tsv', sep='\t')
print(f"Available regions: {len(labels)}")
print(labels[['index', 'name']].head())
```

## Quantitative Results

### Regional Statistics

**Generated by Analysis:**
- Mean field fraction per brain region
- Total blob count per region
- Regional volume measurements  
- Statistical comparisons across subjects

**Example Analysis:**
```python
import nibabel as nib
import pandas as pd

# Load field fraction map and atlas
fieldfrac = nib.load('sub-001_seg-roi82_stain-abeta_fieldfrac.nii').get_fdata()
atlas = nib.load('tpl-ABAv3_seg-roi82_dseg.nii.gz').get_fdata()
labels = pd.read_csv('tpl-ABAv3_seg-roi82_dseg.tsv', sep='\t')

# Calculate regional statistics
results = []
for idx, row in labels.iterrows():
    region_mask = atlas == row['index']
    if region_mask.sum() > 0:
        region_ff = fieldfrac[region_mask]
        results.append({
            'region': row['name'],
            'mean_fieldfrac': region_ff.mean(),
            'std_fieldfrac': region_ff.std(),
            'volume_voxels': region_mask.sum()
        })

regional_stats = pd.DataFrame(results)
print(regional_stats.head())
```

## Quality Control Outputs

### Processing Logs

**Location:**
```
output/spimquant/logs/
├── import/                             # Data import logs
├── templatereg/                        # Registration logs
├── segmentation/                       # Segmentation logs
└── quantification/                     # Analysis logs
```

**Content:**
- Command execution details
- Processing times and resource usage
- Error messages and warnings
- Parameter settings used

### Registration Quality Metrics

**Metrics Available:**
```yaml
Cross-correlation:     # Image similarity measure
  typical_range: 0.6 - 0.9
  interpretation: Higher values indicate better alignment

Mutual_information:    # Information-theoretic measure
  typical_range: 1.0 - 3.0  
  interpretation: Higher values indicate better registration

Jacobian_statistics:   # Deformation field analysis
  mean_expansion: Should be close to 1.0
  std_expansion: Lower values indicate smoother deformation
```

### Visual QC Outputs

**Registration Overlays:**
Generate composite images for visual inspection:
```bash
# Create overlay images
c3d tpl-ABAv3_anat.nii.gz \
    sub-001_space-ABAv3_stain-PI_SPIM.nii \
    -overlay-label-image \
    -o registration_qc.nii.gz
```

## Advanced Output Analysis

### Multi-Subject Analysis

**Aggregation Example:**
```python
import os
import pandas as pd
import nibabel as nib

def extract_regional_data(subject_dir, atlas='roi82'):
    """Extract regional statistics for one subject."""
    fieldfrac_file = f"{subject_dir}/sub-{subject}_seg-{atlas}_stain-abeta_fieldfrac.nii"
    
    if os.path.exists(fieldfrac_file):
        fieldfrac = nib.load(fieldfrac_file).get_fdata()
        atlas_img = nib.load(f'tpl-ABAv3_seg-{atlas}_dseg.nii.gz').get_fdata()
        
        # Calculate regional means
        regional_data = {}
        for region_id in range(1, atlas_img.max() + 1):
            mask = atlas_img == region_id
            if mask.sum() > 0:
                regional_data[f'region_{region_id}'] = fieldfrac[mask].mean()
        
        return regional_data
    return {}

# Process all subjects
all_subjects = []
for subdir in os.listdir('output/spimquant/'):
    if subdir.startswith('sub-'):
        subject_data = extract_regional_data(f'output/spimquant/{subdir}/micr')
        subject_data['subject'] = subdir
        all_subjects.append(subject_data)

# Create group-level dataframe
group_data = pd.DataFrame(all_subjects)
print(f"Group analysis: {len(group_data)} subjects, {len(group_data.columns)-1} regions")
```

### Statistical Analysis

**Region-wise Comparisons:**
```python
from scipy import stats
import matplotlib.pyplot as plt

# Compare groups
group1_data = group_data[group_data['group'] == 'control']['region_5'].dropna()
group2_data = group_data[group_data['group'] == 'treatment']['region_5'].dropna()

# Statistical test
statistic, p_value = stats.ttest_ind(group1_data, group2_data)
print(f"Region 5 comparison: t={statistic:.3f}, p={p_value:.4f}")

# Visualization
plt.figure(figsize=(8, 6))
plt.boxplot([group1_data, group2_data], labels=['Control', 'Treatment'])
plt.ylabel('Field Fraction')
plt.title('Amyloid-beta Field Fraction - Hippocampus')
plt.show()
```

## Data Export and Sharing

### Export to Standard Formats

**CSV Export:**
```python
# Export regional data
regional_summary = group_data.drop('subject', axis=1)
regional_summary.to_csv('spimquant_regional_results.csv', index=False)
```

**NIfTI Export:**
Original NIfTI files are directly compatible with most neuroimaging software:
- FSL
- AFNI  
- SPM
- ANTs
- ITK-SNAP

### Integration with Other Tools

**FSL Integration:**
```bash
# View in FSLEyes
fsleyes tpl-ABAv3_anat.nii.gz \
         sub-001_space-ABAv3_stain-abeta_SPIM.nii \
         -cm hot -dr 0 100
```

**Python Neuroimaging:**
```python
from nilearn import plotting
import matplotlib.pyplot as plt

# Create publication-quality figures
plotting.plot_stat_map(
    'sub-001_seg-roi82_stain-abeta_fieldfrac.nii',
    bg_img='tpl-ABAv3_anat.nii.gz',
    title='Amyloid-beta Field Fraction',
    colorbar=True
)
plt.show()
```

## Troubleshooting Output Issues

### Common Problems

**Missing Output Files:**
- Check processing logs for errors
- Verify input data quality
- Review configuration parameters

**Unexpected Results:**
- Validate registration quality visually
- Check segmentation parameters
- Compare with known standards

**Performance Issues:**
- Monitor disk space usage
- Check memory consumption logs
- Consider reducing resolution levels

### Output Validation

**Automated Checks:**
```bash
# Verify output completeness
find output/spimquant -name "*.nii" | wc -l
find output/spimquant -name "*fieldfrac*" | head -5

# Check file sizes
du -sh output/spimquant/sub-*/micr/*.nii
```

**Visual Validation:**
Always perform visual quality control of key outputs:
1. Registration accuracy
2. Segmentation quality  
3. Template alignment
4. Regional boundary definitions