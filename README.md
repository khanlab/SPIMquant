# SPIMquant
[![Documentation Status](https://readthedocs.org/projects/spimquant/badge/?version=latest)](https://spimquant.readthedocs.io/en/latest/?badge=latest)

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

SPIMquant supports group-level statistical analysis using a formula-based OLS model (via [statsmodels](https://www.statsmodels.org/) / [patsy](https://patsy.readthedocs.io/)), comparing segmentation statistics (e.g., fieldfrac, density, volume) across groups of participants.

1. First, create a `participants.tsv` file in your BIDS directory with participant metadata including group assignment columns (e.g., `treatment`, `genotype`, `sex`, `age`):
   ```tsv
   participant_id	treatment	genotype	sex	age
   sub-01	vehicle	WT	M	12
   sub-02	vehicle	WT	F	13
   sub-03	drug	WT	M	11
   sub-04	drug	WT	F	12
   sub-05	vehicle	KO	M	12
   sub-06	drug	KO	M	11
   ```

2. Run group-level analysis specifying the statistical model and pairwise contrasts:
   ```bash
   pixi run spimquant /path/to/bids/dir /path/to/output/dir group \
     --group-stats-model "metric ~ C(treatment) + age" \
     --group-stats-pairwise treatment \
     --cores all
   ```

   For more complex designs with interaction effects and stratified contrasts:
   ```bash
   pixi run spimquant /path/to/bids/dir /path/to/output/dir group \
     --group-stats-model "metric ~ C(treatment) * C(genotype) * C(sex) + age" \
     --group-stats-pairwise treatment \
     --group-stats-within genotype sex \
     --cores all
   ```

   To restrict the analysis cohort (e.g., exclude subjects not meeting QC criteria):
   ```bash
   pixi run spimquant /path/to/bids/dir /path/to/output/dir group \
     --group-stats-model "metric ~ C(treatment) + age" \
     --group-stats-pairwise treatment \
     --group-stats-where "treatment in ['vehicle', 'drug'] and qc_pass == 1" \
     --cores all
   ```

   Use `--group-stats-label` to name the analysis run so multiple analyses don't overwrite each other (default: `1`):
   ```bash
   pixi run spimquant /path/to/bids/dir /path/to/output/dir group \
     --group-stats-label treatment_analysis \
     --group-stats-model "metric ~ C(treatment) + age" \
     --group-stats-pairwise treatment \
     --cores all
   ```

This will generate (under `<output_dir>/group/<label>/`):
- `*_allsubjects.tsv`: Merged ROI-level table for all subjects (with participant metadata) -- export to your own stats tools
- `*_groupstats.tsv`: Statistical results per region (t-statistic, p-value, Cohen's d, group means) for each pairwise contrast
- `*_groupstats.png`: Heatmap visualizations of statistical results
- `*_groupstats.nii`: 3D volumetric maps of statistical values for visualization in neuroimaging software

# Contributing
 We welcome contributions! Please refer to the [contributing guidelines](CONTRIBUTING.md) for more details on how to contribute.

# License
 This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
