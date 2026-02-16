# Import and Format Conversion

The import workflow handles loading reference data (templates, atlases, masks) and converting SPIM data between OME-Zarr and NIfTI formats.

## Purpose

This workflow module provides the foundation for all downstream processing by:
- Importing template anatomies and atlases from the resources directory
- Converting multi-resolution OME-Zarr data to NIfTI at specific pyramid levels
- Organizing all files according to BIDS conventions
- Supporting both local files and remote URLs for reference data

## Workflow Location

**File**: `spimquant/workflow/rules/import.smk`

## Key Rule Types

### Template Import Rules

These rules import reference templates and associated files:

- **Template Anatomy**: Imports the reference anatomical image (e.g., ABAv3, gubra)
- **Brain Masks**: Imports template brain masks used as registration priors
- **Atlas Segmentations**: Imports parcellation files (dseg) defining brain regions
- **Label Tables**: Imports TSV files mapping region IDs to names and colors

All import rules are marked as `localrules` for efficient execution without cluster submission.

### Format Conversion Rules

- **OME-Zarr to NIfTI**: Extracts specific resolution levels from multi-scale pyramids
- **Label Format Conversion**: Converts TSV label tables to visualization tool formats (e.g., ITK-SNAP)

## Data Types and Formats

### Inputs
- **OME-Zarr** (`.ome.zarr`, `.ome.zarr.zip`): Multi-resolution SPIM data with OME metadata
- **NIfTI** (`.nii.gz`): Template anatomies, masks, and atlases
- **TSV/CSV** (`.tsv`, `.csv`): Label lookup tables

### Outputs
- **NIfTI** (`.nii.gz`): Downsampled SPIM data, imported templates/atlases
- **TSV** (`.tsv`): Standardized label tables
- **ITK-SNAP LUT** (`.itksnap.txt`): Label files for ITK-SNAP visualization

## Resolution Levels

The OME-Zarr format stores data as a multi-resolution pyramid:
- **Level 0**: Highest resolution (e.g., 1.8 μm isotropic)
- **Level 1**: 2x downsampled (e.g., 3.6 μm)
- **Level 2**: 4x downsampled (e.g., 7.2 μm)
- **Level N**: Further downsampling levels

Different processing stages use different resolution levels:
- **Registration**: Typically level 2-3 for computational efficiency
- **Segmentation**: Level 0-1 for detecting small objects
- **Visualization**: Level 0 for high-quality outputs

## Configuration

### Template Configuration

Templates are defined in the configuration file:

```yaml
templates:
  ABAv3:
    anat: "tpl-ABAv3/tpl-ABAv3_res-25um_T2w.nii.gz"
    mask: "tpl-ABAv3/tpl-ABAv3_res-25um_desc-brain_mask.nii.gz"
    atlases:
      AllenCCFv3:
        dseg: "tpl-ABAv3/tpl-ABAv3_res-25um_atlas-AllenCCFv3_dseg.nii.gz"
        tsv: "tpl-ABAv3/tpl-ABAv3_atlas-AllenCCFv3_dseg.tsv"
```

### Orientation Handling

SPIM data may have different orientations than standard NIfTI:

```yaml
orientation: "IPL"  # Specify orientation for proper conversion
```

## Storage Plugins

SPIMquant supports reading data directly from cloud storage:

```python
# S3 storage
inputs["spim"].path  # Automatically handles s3://bucket/path URIs

# Google Cloud Storage
inputs["spim"].path  # Automatically handles gs://bucket/path URIs
```

The `storage()` function wraps paths to enable transparent cloud access.

## BIDS Organization

All imported files follow BIDS conventions:

```
derivatives/
├── tpl-ABAv3/
│   ├── tpl-ABAv3_T2w.nii.gz           # Template anatomy
│   ├── tpl-ABAv3_desc-brain_mask.nii.gz  # Brain mask
│   ├── atlas-AllenCCFv3/
│   │   ├── tpl-ABAv3_atlas-AllenCCFv3_dseg.nii.gz  # Segmentation
│   │   └── tpl-ABAv3_atlas-AllenCCFv3_dseg.tsv     # Label table
```

## Dependencies

- **zarrnii**: OME-Zarr to NIfTI conversion
- **nibabel**: NIfTI file handling
- **pandas**: TSV/CSV table processing
- **s3fs/gcsfs**: Cloud storage access (optional)

## Common Workflows

### Using Different Resolution Levels

```bash
# High-resolution segmentation (level 0)
spimquant /path/to/data /path/to/output participant \
  --segmentation_level 0

# Fast registration at lower resolution (level 2)
spimquant /path/to/data /path/to/output participant \
  --registration_level 2
```

### Adding New Templates

1. Place template files in `spimquant/resources/`
2. Add template configuration to `snakebids.yml`
3. Reference by name in workflow execution

### Remote Data Access

```bash
# Process data directly from S3
spimquant s3://bucket/bids_dataset /path/to/output participant
```

## Related Workflows

- **[Masking](masking.md)**: Uses imported templates for brain extraction
- **[Registration](templatereg.md)**: Registers to imported templates
- **[Segmentation](segmentation.md)**: Uses imported atlases for quantification
