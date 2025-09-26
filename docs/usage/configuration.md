# Configuration

SPIMquant uses YAML configuration files to control processing parameters, template selection, and analysis options. This page describes all available configuration options.

## Configuration File Location

**Default Configuration:**
SPIMquant includes a built-in configuration file at:
```
spimquant/config/snakebids.yml
```

**Custom Configuration:**
You can override settings by providing a custom configuration file:
```bash
spimquant /path/to/data /path/to/output participant --config-file custom_config.yml
```

## Core Parameters

### Analysis Levels
```yaml
analysis_levels:
  - participant      # Process individual subjects

targets_by_analysis_level:
  participant:
    - ''            # Run all rules (empty string = default)
```

### Input Configuration

#### SPIM Data Input
```yaml
pybids_inputs:
  spim:
    filters:
      suffix: 'SPIM'                    # Required BIDS suffix
      extension: 'ome.zarr'             # File extension (.ome.zarr or .ome.zarr.zip)
      sample: brain                     # Sample type filter
    wildcards:
      - subject                         # Subject identifier
      - sample                          # Sample type (e.g., brain)
      - acquisition                     # Acquisition parameters
      - staining                        # Staining protocol
```

#### MRI Data Input (Optional)
```yaml
  T2w:
    filters:
      suffix: 'T2w'                     # T2-weighted MRI
      extension: 'nii.gz'               # NIfTI format
      datatype: 'anat'                  # Anatomical data
    wildcards:
      - subject
      - session
      - acquisition
```

## Processing Parameters

### Template Selection
```yaml
# Primary template for SPIM registration
--template: ABAv3                       # Options: ABAv3, gubra, MBMv3, turone

# MRI template (if using --register_to_mri)
--template_mri: MouseIn                 # Options: MouseIn

# Atlas segmentations to use
--atlas_segs: null                      # null = use all available, or specify list
# Example: --atlas_segs roi22 roi82

# Template cropping for hemisphere-specific analysis
--template_crop: null                   # Options: left, right, null
```

### Registration Parameters
```yaml
# Resolution levels (0 = full resolution, higher = more downsampled)
--registration_level: 5                 # Level for template registration
--segmentation_level: 0                 # Level for final segmentation

# Stains suitable for registration (priority order)
--stains_for_reg:
  - PI                                  # Propidium Iodide
  - YOPRO                              # YoPro-1
  - YoPro                              # Alternative spelling
  - AutoF                              # Autofluorescence
  - autof                              # Lowercase variant

# Stains for quantitative analysis
--stains_for_seg:
  - abeta                              # Amyloid-beta plaques
  - Abeta                              # Alternative spelling
  - BetaAmyloid                        # Full name
  - AlphaSynuclein                     # Alpha-synuclein deposits
  - Iba1                               # Microglial marker
  - ChAT                               # Cholinergic marker
```

### Image Processing
```yaml
# Intensity correction method
--correction_method: n4                 # Options: n4, gaussian

# Segmentation method for pathology detection
--seg_method: otsu+k3i2                # Options: threshold, otsu+k3i2
--seg_threshold: 75                    # Used only with 'threshold' method

# Quality control
--sloppy: false                        # Use low-quality settings for testing
```

### Advanced Options
```yaml
# MRI-based processing
--register_to_mri: false               # Register to subject's MRI first

# Template masking
--template_negative_mask: placeholder   # Path to exclusion mask

# BIDS validation
--skip-bids-validation: false          # Skip input validation
```

## Template Configurations

### ABAv3 (Allen Brain Atlas v3)
```yaml
templates:
  ABAv3:
    anat: '{workflow.basedir}/../resources/ABAv3/P56_Atlas.nii.gz'
    dseg: '{workflow.basedir}/../resources/ABAv3/P56_Annotation.nii.gz'
    lut: '{workflow.basedir}/../resources/ABAv3/labelmapper_ABAv3_to_all.json'
    segs:
      all:                             # Complete atlas
        dseg: 'tpl-ABAv3/tpl-ABAv3_desc-LR_dseg.nii.gz'
        tsv: 'tpl-ABAv3/tpl-ABAv3_desc-LR_dseg.tsv'
      roi22:                           # 22 major regions
        dseg: '{workflow.basedir}/../resources/ABAv3/eed_labels/P56_annotation_22_R_L.nii.gz'
        csv: '{workflow.basedir}/../resources/ABAv3/eed_labels/P56_annotation_22_R_L.csv'
      roi82:                           # 82 brain regions
        dseg: '{workflow.basedir}/../resources/ABAv3/eed_labels/P56_annotation_82_R_L.nii.gz'
        csv: '{workflow.basedir}/../resources/ABAv3/eed_labels/P56_annotation_82_R_L.csv'
      roi198:                          # 198 detailed regions
        dseg: '{workflow.basedir}/../resources/ABAv3/ngo_labels/atlas.nii.gz'
        csv: '{workflow.basedir}/../resources/ABAv3/ngo_labels/atlas_with_names.csv'
```

### Gubra Template
```yaml
  gubra:
    anat: '{workflow.basedir}/../resources/gubra/gubra_template_olf_spacing_reslice.nii.gz'
    dseg: '{workflow.basedir}/../resources/gubra/gubra_ano_olf_spacing_remap_reslice.nii.gz'
    lut: '{workflow.basedir}/../resources/ABAv3/labelmapper_ABAv3_to_all.json'
    segs:
      all:
        dseg: 'tpl-gubra/tpl-gubra_desc-LR_dseg.nii.gz'
        tsv: 'tpl-gubra/tpl-gubra_desc-LR_dseg.tsv'
```

### MBMv3 (Mouse Brain Mapping v3)
```yaml
  MBMv3:
    anat: '{workflow.basedir}/../resources/MBMv3/template_T2w_brain.nii.gz'
    dseg: '{workflow.basedir}/../resources/MBMv3/segmentation_three_types_seg.nii.gz'
    segs:
      paxinos:
        dseg: '{workflow.basedir}/../resources/MBMv3/atlas_MBM_cortex_vPaxinos.nii.gz'
        itksnap: '{workflow.basedir}/../resources/MBMv3/atlas_MBM_cortex_vPaxinos.itksnap.txt'
```

### Turone Template
```yaml
  turone:
    anat: '{workflow.basedir}/../resources/Turone_Mouse_Brain_Template/Turone_Mouse_Brain_Template/TMBTA_Brain_Template.nii.gz'
    dseg: '{workflow.basedir}/../resources/Turone_Mouse_Brain_Template/Turone_Mouse_Brain_Template/TMBTA_tissue_dseg.nii.gz'
    segs:
      all:
        dseg: '{workflow.basedir}/../resources/Turone_Mouse_Brain_Template/Turone_Mouse_Brain_Atlas/TMBTA_Brain_Atlas.nii.gz'
        itksnap: '{workflow.basedir}/../resources/Turone_Mouse_Brain_Template/Turone_Mouse_Brain_Atlas/TMBTA_ItK_Label_File.txt'
```

## Processing Settings

### OME-Zarr Configuration
```yaml
ome_zarr:
  max_downsampling_layers: 4           # Resolution levels: 0,1,2,3,4
  rechunk_size:                        # Chunk size for processing
    - 1                                # Z dimension
    - 1024                             # Y dimension  
    - 1024                             # X dimension
  scaling_method: 'local_mean'         # Options: nearest, gaussian, local_mean, zoom
```

### Registration Settings
```yaml
templatereg:
  desc: N4brain                        # Processing description
  zooms:
    - 40                               # Zoom levels for registration
```

### Masking Parameters
```yaml
masking:
  gmm_k: 9                            # Gaussian mixture components
  gmm_bg_class: 1                     # Background class
  pre_atropos_downsampling: '50%'     # Downsampling for initial segmentation
```

### Blob Detection
```yaml
blobdetect:
  level: 5                            # Downsampling level for detection
  dseg_level: 5                       # Template downsampling for labeling
  dseg_template: ABAv3                # Template for blob labeling
```

## Example Configurations

### High-Resolution Processing
```yaml
# For high-resolution detailed analysis
registration_level: 3                 # Less downsampling
segmentation_level: 0                 # Full resolution
correction_method: n4                 # Best quality correction
seg_method: otsu+k3i2                # Adaptive thresholding
```

### Fast/Testing Configuration
```yaml
# For quick testing or low-resource systems
registration_level: 6                 # More downsampling
segmentation_level: 3                 # Moderate resolution
sloppy: true                         # Enable fast mode
correction_method: gaussian          # Faster correction
```

### Multi-Stain Pathology Analysis
```yaml
# Optimized for pathology studies
stains_for_seg:
  - abeta                            # Amyloid plaques
  - AlphaSynuclein                   # Synuclein deposits
  - Iba1                             # Microglia
atlas_segs:                          # Multiple resolution levels
  - roi22                            # Major regions
  - roi82                            # Detailed regions
template: ABAv3                      # Standard reference
```

### Cloud Processing Setup
```yaml
# For cloud-based workflows
work: /tmp/spimquant_work            # Fast temporary storage
root: /output/spimquant              # Cloud storage output
use_apptainer: true                  # Container isolation
cores: 16                            # Parallel processing
resources:
  mem_mb: 32000                      # Memory per job
```

## Parameter Validation

**Required Parameters:**
- `bids_dir_or_uri`: Input BIDS dataset path or cloud URI
- Analysis level must be 'participant'
- At least one template must be available

**Parameter Interactions:**
- `--register_to_mri` requires T2w input data
- `--template_crop` only works with symmetric templates
- `--seg_method threshold` requires `--seg_threshold`
- `--atlas_segs` must match available segmentations for chosen template

**Resource Considerations:**
- Higher resolution levels (lower numbers) require more memory
- Multiple stains increase processing time linearly
- Template registration is the most memory-intensive step
- Blob detection scales with segmentation complexity