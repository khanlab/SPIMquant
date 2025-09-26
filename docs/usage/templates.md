# Templates and Data

SPIMquant supports multiple brain templates and data formats for flexible analysis of SPIM microscopy data. This section describes available templates, data format requirements, and how to work with different input types.

## Supported Templates

### ABAv3 (Allen Brain Atlas v3) - Default

**Description:**
The Allen Brain Atlas version 3 is the default template, providing comprehensive coverage of mouse brain anatomy with multiple resolution levels.

**Features:**
- **Age**: P56 (adult mouse)
- **Resolution**: 25μm isotropic
- **Coverage**: Full brain including olfactory regions
- **Coordinate System**: RAS+ orientation
- **Atlas Levels**: 3 different parcellation resolutions

**Atlas Resolutions:**
```yaml
roi22:   22 major brain regions (cortex, hippocampus, thalamus, etc.)
roi82:   82 brain regions (detailed cortical and subcortical areas)
roi198:  198 fine-grained regions (maximum anatomical detail)
```

**Use Cases:**
- Standard mouse brain studies
- Cross-study comparisons
- Publications requiring ABA compatibility
- Multi-resolution analysis needs

### Gubra Template

**Description:**
High-resolution template from Gubra-ApS, optimized for lightsheet microscopy registration.

**Features:**
- **Age**: Adult mouse
- **Resolution**: 25μm isotropic  
- **Coverage**: Full brain with olfactory bulb
- **Optimization**: Specifically designed for SPIM data
- **Atlas**: Uses ABA-compatible labeling

**Advantages:**
- Better contrast for lightsheet data
- Optimized for autofluorescence registration
- Validated for SPIM applications

**Use Cases:**
- Lightsheet-specific studies
- Autofluorescence-based registration
- High-contrast template needs

### MBMv3 (Mouse Brain Mapping v3)

**Description:**
High-field MRI-based template with detailed cortical parcellation.

**Features:**
- **Age**: Adult mouse
- **Modality**: High-resolution MRI
- **Resolution**: Variable (optimized per region)
- **Atlas**: Paxinos-based parcellation
- **Specialization**: Cortical area definitions

**Use Cases:**
- Cortical analysis studies
- MRI-SPIM correlation
- Paxinos atlas compatibility

### Turone Template

**Description:**
Template and atlas from Turone laboratory with tissue-type segmentation.

**Features:**
- **Coverage**: Full mouse brain
- **Atlas**: Anatomical region definitions
- **Labels**: ITK-SNAP compatible format
- **Segmentation**: Tissue type classification

**Use Cases:**
- Tissue-type specific analysis
- ITK-SNAP visualization
- Alternative parcellation scheme

## Data Format Requirements

### Input Data Formats

#### OME-Zarr (Recommended)

**Directory Structure:**
```
sub-001_sample-brain_stain-PI_SPIM.ome.zarr/
├── .zgroup
├── .zattrs
├── 0/           # Full resolution
├── 1/           # 2x downsampled
├── 2/           # 4x downsampled
└── ...
```

**Zip Archives:**
```bash
# Also supported
sub-001_sample-brain_stain-PI_SPIM.ome.zarr.zip
```

**Advantages:**
- Multi-resolution pyramids
- Chunked storage for efficient access
- Metadata preservation
- Cloud-optimized format

**Command Line Usage:**
```bash
# For zip archives
spimquant /path/to/data /output participant --filter-spim extension='ome.zarr.zip'

# For directories
spimquant /path/to/data /output participant --filter-spim extension='ome.zarr'
```

#### BIDS Dataset Structure

**Required Structure:**
```
dataset/
├── dataset_description.json
├── participants.tsv
├── sub-001/
│   └── micr/
│       ├── sub-001_sample-brain_stain-PI_SPIM.ome.zarr/
│       ├── sub-001_sample-brain_stain-abeta_SPIM.ome.zarr/
│       └── sub-001_sample-brain_stain-Iba1_SPIM.ome.zarr/
└── sub-002/
    └── micr/
        └── ...
```

**BIDS Entities:**
- `sub-<label>`: Subject identifier
- `sample-<label>`: Sample type (typically 'brain')
- `stain-<label>`: Fluorescent stain/marker
- `acq-<label>`: Acquisition parameters (optional)

### Stain Naming Conventions

#### Registration Stains
Suitable for template registration (nuclear/structural):
```yaml
Recommended names:
  - PI          # Propidium Iodide
  - YOPRO       # YoPro-1 nucleic acid stain
  - DAPI        # DAPI nuclear stain  
  - AutoF       # Autofluorescence
  - autof       # Autofluorescence (lowercase)
  - Hoechst     # Hoechst nuclear stain
```

#### Analysis Stains  
For pathology quantification:
```yaml
Pathology markers:
  - abeta           # Amyloid-beta plaques
  - Abeta           # Alternative spelling
  - BetaAmyloid     # Full name
  - AlphaSynuclein  # Alpha-synuclein deposits
  - pTau            # Phosphorylated tau
  - Tau             # Total tau

Cellular markers:
  - Iba1            # Microglia/macrophages
  - GFAP            # Astrocytes
  - NeuN            # Neurons
  - ChAT            # Cholinergic neurons
  - TH              # Tyrosine hydroxylase
  - Olig2           # Oligodendrocytes
```

## Cloud Storage Support

### Supported Backends

#### Amazon S3
```bash
# S3 URI format
spimquant s3://bucket-name/bids-dataset /output participant

# With credentials
export AWS_ACCESS_KEY_ID=your_key
export AWS_SECRET_ACCESS_KEY=your_secret
spimquant s3://bucket-name/bids-dataset /output participant
```

#### Google Cloud Storage  
```bash
# GCS URI format
spimquant gs://bucket-name/bids-dataset /output participant

# With service account
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/service-account.json
spimquant gs://bucket-name/bids-dataset /output participant
```

#### Azure Blob Storage
```bash
# Azure URI format  
spimquant az://container-name/bids-dataset /output participant

# With credentials
export AZURE_STORAGE_ACCOUNT=your_account
export AZURE_STORAGE_KEY=your_key
```

### Cloud Configuration

**Storage Plugins:**
SPIMquant automatically uses appropriate storage plugins:
- `snakemake-storage-plugin-s3` for S3
- `snakemake-storage-plugin-gcs` for Google Cloud
- `snakemake-storage-plugin-azure` for Azure

**Performance Optimization:**
```yaml
# Increase concurrent downloads
storage:
  s3:
    max_requests_per_second: 100
  gcs:  
    cache_size: 1000
```

## Working with Templates

### Template Selection

**Command Line:**
```bash
# Choose template
spimquant /data /output participant --template ABAv3

# Choose specific atlas resolution
spimquant /data /output participant --template ABAv3 --atlas_segs roi82

# Multiple atlas levels
spimquant /data /output participant --atlas_segs roi22 roi82 roi198
```

**Configuration File:**
```yaml
template: ABAv3
atlas_segs:
  - roi82
  - roi198
```

### Template Cropping

**Hemisphere Analysis:**
```bash
# Left hemisphere only
spimquant /data /output participant --template_crop left

# Right hemisphere only  
spimquant /data /output participant --template_crop right
```

**Use Cases:**
- Unilateral lesion studies
- Hemisphere-specific analysis
- Reduced computational requirements

### Custom Templates

**Adding New Templates:**
1. **Prepare Template Files:**
   ```
   resources/
   └── tpl-MyTemplate/
       ├── tpl-MyTemplate_anat.nii.gz      # T2-like anatomy
       ├── tpl-MyTemplate_dseg.nii.gz      # Segmentation
       └── tpl-MyTemplate_dseg.tsv         # Region labels
   ```

2. **Update Configuration:**
   ```yaml
   templates:
     MyTemplate:
       anat: '{workflow.basedir}/../resources/tpl-MyTemplate/tpl-MyTemplate_anat.nii.gz'
       dseg: '{workflow.basedir}/../resources/tpl-MyTemplate/tpl-MyTemplate_dseg.nii.gz'
       segs:
         all:
           dseg: 'tpl-MyTemplate/tpl-MyTemplate_dseg.nii.gz'
           tsv: 'tpl-MyTemplate/tpl-MyTemplate_dseg.tsv'
   ```

3. **Label File Format (TSV):**
   ```tsv
   index	name	abbreviation	color
   1	Cortex	Ctx	255,0,0
   2	Hippocampus	Hip	0,255,0
   3	Thalamus	Thal	0,0,255
   ```

## Quality Control

### Template Registration QC

**Visual Inspection:**
```bash
# Load registered image and template
itksnap \
  -g output/tpl-ABAv3/tpl-ABAv3_anat.nii.gz \
  -o output/sub-001/micr/sub-001_space-ABAv3_stain-PI_SPIM.nii
```

**Automated Metrics:**
- Cross-correlation coefficient
- Mutual information
- Dice overlap with brain mask
- Jacobian determinant statistics

### Data Format Validation

**OME-Zarr Inspection:**
```python
import zarr
from zarrnii import ZarrNii

# Load and inspect
znii = ZarrNii.from_ome_zarr('path/to/file.ome.zarr')
print(f"Shape: {znii.shape}")
print(f"Channels: {znii.list_channels()}")
print(f"Resolution levels: {znii.n_levels}")
```

**BIDS Validation:**
```bash
# Validate BIDS structure
bids-validator /path/to/bids/dataset

# Skip validation in SPIMquant
spimquant /data /output participant --skip-bids-validation
```

## Troubleshooting

### Common Issues

**Registration Failures:**
- Check stain selection for registration
- Verify template orientation
- Consider template cropping for partial data

**Memory Issues:**
- Increase registration_level (more downsampling)
- Reduce number of parallel jobs
- Use smaller chunk sizes

**Format Errors:**
- Validate OME-Zarr structure
- Check BIDS naming conventions
- Verify file permissions

**Template Mismatches:**
- Ensure proper species/age matching
- Check coordinate system orientation
- Validate resolution compatibility