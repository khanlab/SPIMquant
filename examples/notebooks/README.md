# SPIMquant Example Notebooks

This directory contains Jupyter notebooks demonstrating various aspects of SPIMquant processing and analysis.

## Available Notebooks

### Colocalization Analysis

**`colocalization_visualization.ipynb`** - Interactive tutorial demonstrating colocalization analysis
- Creates toy regionprops data for two stains (Abeta and Iba1)
- Visualizes spatial distribution with custom colors (magenta for Abeta, yellow with black outline for Iba1)
- Demonstrates KDTree-based spatial indexing for efficient nearest neighbor search
- Shows distance calculations and overlap ratio computation
- Provides multiple visualization types:
  - 2D scatter plots of object positions
  - Nearest neighbor search with distance annotations
  - Overlap ratio calculation examples
  - Colocalization result maps with spatial links
  - Statistical distributions (distance, overlap ratio)
- Includes detailed explanations of the colocalization algorithm used in SPIMquant

**Use case:** Understanding the colocalization analysis approach, creating explanatory figures for presentations

### Image Processing

**`vis_blob_detection.ipynb`** - Visualize blob detection results

**`otsu.ipynb`** - Otsu thresholding demonstration

**`quantify_BetaAmyloid.ipynb`** - Beta-amyloid quantification example

### Maximum Intensity Projection (MIP)

**`MIP.ipynb`** - Basic maximum intensity projection

**`MIP_coiled.ipynb`** - MIP using Coiled for cloud execution

**`MIP_rot3d.ipynb`** - 3D rotation MIP

### Visualization and Analysis

**`example_zarr_process.ipynb`** - Processing OME-Zarr files

**`hist_summary.ipynb`** - Histogram summary statistics

## Running the Notebooks

### Using pixi (recommended)

```bash
# Install with development environment (includes jupyterlab)
pixi install --environment dev

# Start JupyterLab
pixi run -e dev jupyter lab
```

### Using pip

```bash
# Install dependencies
pip install numpy pandas scipy matplotlib jupyterlab

# Start JupyterLab
jupyter lab
```

## Dependencies

Most notebooks require:
- numpy
- pandas
- scipy
- matplotlib

Some notebooks may require additional dependencies:
- zarr, zarrnii (for OME-Zarr processing)
- dask (for parallel processing)
- napari (for visualization)

## Creating New Notebooks

When creating new example notebooks:
1. Add clear markdown documentation explaining the purpose
2. Include a brief description in this README
3. Test that the notebook runs from a fresh kernel
4. Consider using toy/synthetic data for pedagogical examples
5. Use appropriate figure sizes and DPI for presentations
