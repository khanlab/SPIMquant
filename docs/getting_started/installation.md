## SPIMquant

SPIMquant is a BIDS App for processing SPIM (lightsheet) microscopy datasets, performing registration to a template, and quantifying microscopic features from the SPIM data.

Hardware requirements: If run locally, make sure you have sufficient memory 
(at least quite a bit more than 16G of memory in total), as the `greedy` diffeormorphic registration we rely on can consume a significant amount of memory during the template registration process.

Software requirements: A linux machine with pixi installed is 
recommended. Pixi will manage all dependencies including Python packages and external tools
like `itk-snap` (`c3d` comes with itk-snap), `greedy` command line tool, and python environment
according to pyproject.toml.

(For developer of SPIMquant: To update dataset, replace the link below)
A sample dataset can be downloaded from [here](https://drive.google.com/file/d/1-eVG_1VREKCE8auI81jyW4onyUcpzXI7/view?usp=sharing) 
which contains:
- An ome-zarr image of down-sampled mousebrain scan in a BIDS dataset
- Reference template files

## Installation

Install pixi (if not already installed):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Clone the repository and install dependencies:

```bash
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant
pixi install
```

Before running the app, you need to specify a config file to use. "SPIMquant/examples/snakebids_template.yml" 
provides a starting point for specifying a config. If you are using the example dataset provided in the 
above section, then a config file is also included in the zip file.

To specify the config, copy the config file into SPIMquant/spimquant/config/snakebids.yml, and change the 
properties in the config file to ensure paths to the directory are properly set.

## Running the app

Do a dry-run first (`-n`) and simply print (`-p`) what would be run:

```bash
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant -np
```

Run the app, using all cores::

```bash
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant --cores all
```


## Generating a report

After your processing is complete, you can use snakemake's `--report` feature to generate
an HTML report. This report will include a graph of all the jobs run, with clickable nodes
to inspect the shell command or python code used in each job, along with the config files and
run times for each job. 

To generate a report, run:

```bash
pixi run spimquant /path/to/bids/dir /path/to/output/dir participant --report
```


