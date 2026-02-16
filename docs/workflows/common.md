# Common Utilities

The common utilities module provides shared helper functions used across all SPIMquant workflows.

## Purpose

This module centralizes common functionality to:
- Ensure consistent file path construction
- Abstract BIDS naming conventions
- Simplify resource file access
- Provide configuration accessors
- Validate data consistency across subjects

These utilities are imported by other workflow modules and enable cleaner, more maintainable code.

## Workflow Location

**File**: `spimquant/workflow/rules/common.smk`

**Note**: This file contains a module-level docstring describing all utility functions.

## Utility Functions

### `bids_tpl()`

**Purpose**: Construct BIDS paths for template-specific files.

**Signature**:
```python
def bids_tpl(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
```

**Usage**:
```python
# Get template anatomy path
template_path = bids_tpl(
    root=root,
    template="ABAv3",
    suffix="anat.nii.gz"
)
# Returns: "{root}/tpl-ABAv3/tpl-ABAv3_anat.nii.gz"

# Get atlas segmentation path
dseg_path = bids_tpl(
    root=root,
    template="ABAv3",
    seg="AllenCCFv3",
    suffix="dseg.nii.gz"
)
# Returns: "{root}/tpl-ABAv3/atlas-AllenCCFv3/tpl-ABAv3_atlas-AllenCCFv3_dseg.nii.gz"
```

### `resources_path()`

**Purpose**: Resolve paths to resource files, supporting both local and remote URLs.

**Signature**:
```python
def resources_path(path):
    """Get path relative to the resources folder"""
```

**Usage**:
```python
# Local file
local_template = resources_path("tpl-ABAv3/template.nii.gz")
# Returns: "{workflow.basedir}/../resources/tpl-ABAv3/template.nii.gz"

# Remote URL (passed through unchanged)
remote_template = resources_path("https://example.com/template.nii.gz")
# Returns: "https://example.com/template.nii.gz"
```

**Behavior**:
- Local paths: Resolved relative to `spimquant/resources/`
- HTTP/HTTPS URLs: Returned unchanged
- Enables flexible resource storage (local or cloud)

### `get_template_path()`

**Purpose**: Get template path with optional hemisphere cropping.

**Signature**:
```python
def get_template_path(root, template, template_crop=None):
    """Get template path, optionally cropped based on hemisphere"""
```

**Usage**:
```python
# Full template
full_path = get_template_path(root, "ABAv3", template_crop=None)

# Left hemisphere only
left_path = get_template_path(root, "ABAv3", template_crop="left")
# Uses: tpl-ABAv3_desc-leftcrop_anat.nii.gz

# Right hemisphere only
right_path = get_template_path(root, "ABAv3", template_crop="right")
# Uses: tpl-ABAv3_desc-rightcrop_anat.nii.gz
```

### `get_template_for_reg()`

**Purpose**: Get the appropriate template for registration based on configuration.

**Signature**:
```python
def get_template_for_reg(wildcards):
    """Get the appropriate template file for registration, cropped if specified"""
```

**Usage**:
Used in rule inputs to dynamically select template:
```python
rule affine_reg:
    input:
        template=get_template_for_reg,
        subject=...,
```

**Behavior**:
- Checks `config["template_crop"]`
- Returns cropped or full template as appropriate
- Ensures consistent template usage across registration rules

### `get_stains_all_subjects()`

**Purpose**: Validate that all subjects have consistent channel/stain names.

**Signature**:
```python
def get_stains_all_subjects():
    """Get stain list, ensuring consistency across subjects"""
```

**Usage**:
```python
# Called during workflow initialization
stains = get_stains_all_subjects()
# Returns: ["PI", "Abeta", "Iba1"]

# Raises ValueError if subjects have different stains
```

**Behavior**:
- Reads channel names from OME-Zarr metadata
- Compares across all subjects
- Raises informative error if inconsistent
- Ensures workflow can safely expand across stains

## Configuration Access Patterns

### Safe Dictionary Access

The configuration is accessed using safe patterns:

```python
# With default value
value = config.get("key", default_value)

# Nested with fallback
value = config.get("section", {}).get("key", default_value)

# Required value (will error if missing)
value = config["required_key"]
```

### Common Configuration Keys

```python
# Template selection
template = config["template"]
template_crop = config.get("template_crop", None)

# Resolution levels
registration_level = config["registration_level"]
segmentation_level = config["segmentation_level"]

# Stain selection
stain_for_reg = config["templatereg"]["stain"]
stains_for_seg = config["stains_for_seg"]

# Processing parameters
sloppy = config.get("sloppy", False)
orientation = config["orientation"]
```

## BIDS Conventions

### Entity Order

SPIMquant follows BIDS entity ordering:
1. Subject: `sub-`
2. Session: `ses-`
3. Sample: `sample-`
4. Acquisition: `acq-`
5. Space: `space-`
6. Template: `tpl-`
7. Atlas/Segmentation: `atlas-`, `seg-`
8. Description: `desc-`
9. Suffix: Last component before extension

### File Naming Examples

```python
# Subject-specific SPIM data
"sub-001_sample-brain_stain-PI_level-2_desc-N4brain_SPIM.nii.gz"

# Template anatomy
"tpl-ABAv3_res-25um_T2w.nii.gz"

# Atlas segmentation
"tpl-ABAv3_atlas-AllenCCFv3_dseg.nii.gz"

# Transformed image
"sub-001_space-ABAv3_stain-PI_desc-deform_SPIM.nii.gz"

# Statistics
"sub-001_seg-AllenCCFv3_from-ABAv3_desc-otsu+k3i2+cleaned_mergedsegstats.tsv"
```

## Directory Structure

### Standard Layout

```
derivatives/
├── sub-{subject}/
│   ├── micr/                    # Microscopy derivatives
│   │   ├── *_SPIM.nii.gz        # Processed images
│   │   ├── *_mask.nii.gz        # Masks
│   │   ├── *_regionprops.parquet
│   │   └── *_segstats.tsv
│   ├── warps/                   # Transformations
│   │   ├── *_xfm.txt
│   │   └── *_warp.nii.gz
│   └── anat/                    # MRI derivatives (if applicable)
├── tpl-{template}/              # Template-specific files
│   ├── tpl-{template}_*.nii.gz
│   └── atlas-{atlas}/
│       └── tpl-{template}_atlas-{atlas}_*.nii.gz
├── group/                       # Group-level results
│   ├── *_groupstats.tsv
│   └── *_groupavg.nii.gz
└── logs/                        # Processing logs
```

## Error Handling

### Common Validation Patterns

```python
# Check for file existence
if not Path(input_file).exists():
    raise FileNotFoundError(f"Input file not found: {input_file}")

# Validate configuration
if "required_key" not in config:
    raise ValueError("Configuration must include 'required_key'")

# Check stain consistency
stain_sets = [get_stains(zarr) for zarr in inputs["spim"]]
if not all(s == stain_sets[0] for s in stain_sets):
    raise ValueError(f"Inconsistent stains across subjects: {stain_sets}")
```

## Usage in Workflows

### Typical Import Pattern

```python
# At top of workflow .smk file
from pathlib import Path
from snakebids import bids

# Common utilities are automatically available
# from common.smk include:
# - bids_tpl()
# - resources_path()
# - get_template_for_reg()
# etc.
```

### Rule Example

```python
rule my_rule:
    input:
        # Use utility functions in input/output
        template=get_template_for_reg,
        dseg=bids_tpl(
            root=root,
            template="{template}",
            seg="{seg}",
            suffix="dseg.nii.gz"
        ),
    output:
        result=bids(
            root=root,
            datatype="micr",
            **inputs["spim"].wildcards,
            suffix="result.nii.gz"
        ),
```

## Best Practices

### Path Construction

**Do**:
```python
# Use bids() and bids_tpl() for all paths
path = bids_tpl(root=root, template="ABAv3", suffix="anat.nii.gz")
```

**Don't**:
```python
# Avoid manual path construction
path = f"{root}/tpl-ABAv3/tpl-ABAv3_anat.nii.gz"  # Fragile!
```

### Configuration Access

**Do**:
```python
# Use .get() with defaults
value = config.get("optional_key", default_value)
```

**Don't**:
```python
# Avoid direct access to optional keys
value = config["optional_key"]  # Raises KeyError if missing
```

### Resource Access

**Do**:
```python
# Use resources_path() for all resources
path = resources_path(config["template"]["anat"])
```

**Don't**:
```python
# Avoid hardcoded paths
path = "../resources/template.nii.gz"  # Not portable!
```

## Testing and Validation

### Verify Path Construction

```python
# Test bids_tpl output
test_path = bids_tpl(root="test", template="ABAv3", suffix="anat.nii.gz")
assert "tpl-ABAv3" in test_path
assert test_path.endswith("anat.nii.gz")
```

### Validate Stain Consistency

```bash
# Dry-run to check for stain errors
spimquant /path/to/data /path/to/output participant -n
```

If subjects have inconsistent stains, workflow will fail early with clear error message.

## Related Workflows

All workflows depend on these utilities:
- **[Import](import.md)**: Uses `bids_tpl()` and `resources_path()`
- **[Masking](masking.md)**: Uses `get_template_for_reg()`
- **[Registration](templatereg.md)**: Uses all utility functions extensively
- **[Segmentation](segmentation.md)**: Uses BIDS path construction
- **[Group Statistics](groupstats.md)**: Uses path utilities for aggregation

## Extension Points

### Adding New Utilities

To add new utility functions:

1. Add to `common.smk`:
```python
def my_new_utility(arg1, arg2):
    """Description of utility"""
    # implementation
    return result
```

2. Document in this file

3. Use in other workflow modules

### Custom BIDS Entities

To support new BIDS entities:

1. Update `bids()` calls to include new entity:
```python
bids(
    root=root,
    custom_entity="{value}",
    suffix="file.nii.gz"
)
```

2. Add to documentation
3. Ensure entity appears in correct order per BIDS specification
