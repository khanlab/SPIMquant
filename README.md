# SPIMquant

SPIMquant is a Snakebids app for quantitative analysis of SPIM (lightsheet) brain data. It performs automated nonlinear template registration and quantification of pathology from SPIM microscopy datasets.

Features include:
 - Deformable registration to a template 
 - Atlas-based quantification of pathology 
 - High-resolution Imaris dataset creation from atlas region bounding boxes
 - Group-level statistical analysis with contrast comparisons
 - Coarse-grained and fine-grained parallelization using Snakemake and Dask
 - Support for reading BIDS datasets directly from cloud-based object storage
 - Support for simple and scalable cloud-based processing with Coiled


# Installation

## Hardware Requirements
 - Processing lightsheet microscopy data is computationally-demanding, and you will need sufficient (and ideally fast and local) 
 disk space. The more cores you have access to, the faster the code will run, but you will also need sufficient memory (e.g. 2-4 GB per core) as well.

## Software Requirements
 - A Linux machine with pixi installed. Pixi will manage all Python dependencies and non-python dependencies (c3d, greedy, ANTS) through conda environments.

## Steps
 1. Install pixi (if not already installed):
    ```bash
    curl -fsSL https://pixi.sh/install.sh | bash
    ```
    
 2. Clone the repository and install dependencies:
    ```bash
    git clone https://github.com/khanlab/SPIMquant.git
    cd SPIMquant
    pixi install
    ```
    
## Development
For development work, use the development environment which includes additional tools like formatters and linters:

```bash
pixi install --environment dev
```

Quality checks can be run using:
```bash
pixi run quality_check  # Check code formatting
pixi run quality_fix    # Fix code formatting
```   

# Usage

SPIMquant is a BIDS App, so you need a BIDS dataset containing SPIM (or lightsheet microscopy) data to use it. The [SPIMprep](https://github.com/khanlab/SPIMprep)
workflow is the recommended tool to produce a BIDS dataset from your raw or minimally-preprocessed microscopy data.

 1. Perform a dry run:
    ```bash
    pixi run spimquant /path/to/bids/dir /path/to/output/dir participant -np
    ```
 2. Run the app using all cores:
    ```bash
    pixi run spimquant /path/to/bids/dir /path/to/output/dir participant --cores all
    ```

If your input BIDS dataset stores data in zarr zipstores (e.g. SPIM files ending in `*_SPIM.ome.zarr.zip`), then you should use the following option:
```
--filter-spim extension='ome.zarr.zip'
```

## Group-Level Statistical Analysis

SPIMquant supports group-level statistical analysis to compare segmentation statistics (e.g., fieldfrac, density, volume) across groups of participants. 

1. First, create a `participants.tsv` file in your BIDS directory with participant metadata including a column for group assignments (e.g., 'treatment', 'genotype'):
   ```tsv
   participant_id	treatment	age	sex
   sub-01	control	12	M
   sub-02	control	13	F
   sub-03	drug	11	M
   sub-04	drug	12	F
   ```

2. Run group-level analysis specifying the contrast column and values:
   ```bash
   pixi run spimquant /path/to/bids/dir /path/to/output/dir group \
     --contrast_column treatment \
     --contrast_values control drug \
     --cores all
   ```

This will generate:
- `*_groupstats.tsv`: Statistical test results (t-statistics, p-values, effect sizes) for each brain region
- `*_groupstats.png`: Heatmap visualizations of statistical results
- `*_groupstats.nii`: 3D volumetric maps of statistical values for visualization in neuroimaging software

# Contributing
 We welcome contributions! Please refer to the [contributing guidelines](CONTRIBUTING.md) for more details on how to contribute.

# License
 This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
