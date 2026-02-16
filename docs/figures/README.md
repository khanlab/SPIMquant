# SPIMquant Workflow DAG Diagrams

This directory contains Mermaid diagrams representing the SPIMquant workflow directed acyclic graph (DAG).

## Overview

The workflow is broken down into **11 functional stages**, each with its own diagram:

1. **01_import** - Import and setup templates, masks, and atlases
2. **02_preprocessing** - Image preprocessing and downsampling
3. **03_masking** - Brain masking using atropos segmentation
4. **04_correction** - Intensity correction and normalization
5. **05_registration** - Template registration (affine and deformable)
6. **06_transform** - Apply transformations to images and atlases
7. **07_segmentation** - Segmentation of pathology (threshold, multi-otsu)
8. **08_quantification** - Region properties and quantification
9. **09_statistics** - Statistical analysis and atlas mapping
10. **10_qc** - Quality control and reporting
11. **11_patches** - Extract image patches for analysis

## Files

- `rulegraph_full.mermaid` - Complete workflow rulegraph showing all rules and dependencies
- `dag_01_import.mermaid` - Import stage diagram
- `dag_02_preprocessing.mermaid` - Preprocessing stage diagram
- `dag_03_masking.mermaid` - Masking stage diagram
- `dag_04_correction.mermaid` - Correction stage diagram
- `dag_05_registration.mermaid` - Registration stage diagram
- `dag_06_transform.mermaid` - Transform stage diagram
- `dag_07_segmentation.mermaid` - Segmentation stage diagram
- `dag_08_quantification.mermaid` - Quantification stage diagram
- `dag_09_statistics.mermaid` - Statistics stage diagram
- `dag_10_qc.mermaid` - QC stage diagram
- `dag_11_patches.mermaid` - Patches stage diagram

## Viewing Diagrams

### Online Viewers

You can view these Mermaid diagrams using:
- [Mermaid Live Editor](https://mermaid.live/) - Paste the diagram content
- [GitHub's native Mermaid support](https://github.blog/2022-02-14-include-diagrams-markdown-files-mermaid/) - View directly in GitHub markdown

### Local Rendering

To render diagrams to SVG/PNG locally:

```bash
# Install mermaid-cli
npm install -g @mermaid-js/mermaid-cli

# Render a single diagram
mmdc -i dag_05_registration.mermaid -o dag_05_registration.svg

# Render all diagrams
for f in *.mermaid; do mmdc -i "$f" -o "${f%.mermaid}.svg"; done
```

## Regenerating Diagrams

To regenerate all diagrams from the current workflow:

```bash
cd docs/scripts
python3 generate_dag_diagrams.py
```

This will:
1. Run SPIMquant on the `tests/bids_ds` dataset to generate the complete DAG
2. Parse the DAG and classify rules into functional stages
3. Create individual Mermaid diagrams for each stage
4. Optionally render to SVG if mermaid-cli is installed

### Advanced Options

```bash
# Use a different test dataset
python3 generate_dag_diagrams.py --bids-dir /path/to/bids/dataset

# Use a different template
python3 generate_dag_diagrams.py --template gubra

# Specify custom output directory
python3 generate_dag_diagrams.py --output-dir /tmp/custom_output
```

## Diagram Legend

In the stage-specific diagrams:

- **Bold green nodes** (thick border) - Primary rules for this stage
- **Light blue nodes** (thin border) - Dependency rules from other stages
- **Arrows** - Dependencies (from source to target)

## CI Integration

The diagram generation script can be integrated into CI/CD pipelines:

```yaml
# Example GitHub Actions workflow
- name: Generate DAG diagrams
  run: |
    python3 docs/scripts/generate_dag_diagrams.py
    
# Optional: Render to SVG
- name: Install mermaid-cli
  run: npm install -g @mermaid-js/mermaid-cli
  
- name: Render diagrams
  run: |
    cd docs/figures
    for f in *.mermaid; do 
      mmdc -i "$f" -o "${f%.mermaid}.svg" -t neutral
    done
```

## Workflow Stages Details

### 01_import
Handles importing template anatomical images, brain masks, atlas segmentations, and label lookup tables. Sets up the reference space for registration.

### 02_preprocessing
Converts OME-Zarr multiscale images to NIfTI format at specified downsampling levels. Prepares data for registration and analysis.

### 03_masking
Creates brain masks using Atropos segmentation with Gaussian mixture models. Applies template masks to subject space as priors.

### 04_correction
Applies N4 bias field correction to reduce intensity non-uniformities in the images. Masks corrected images to brain tissue.

### 05_registration
Performs multi-stage registration: initialization, affine, and deformable registration to align subject images to template space.

### 06_transform
Applies computed transformations to warp images, segmentations, and region properties between subject and template spaces.

### 07_segmentation
Segments pathology features using thresholding and multi-Otsu methods. Computes field fractions and density maps.

### 08_quantification
Extracts region properties from segmented objects, filters based on size criteria, and maps results to atlas regions.

### 09_statistics
Aggregates statistics across atlas regions, merges results from multiple stains, and creates quantitative feature maps.

### 10_qc
Generates quality control reports with registration overlays and segmentation visualizations for manual inspection.

### 11_patches
Extracts 3D image patches from specific atlas regions for downstream machine learning or detailed analysis.
