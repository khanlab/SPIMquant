## The `spimquant` Interface

The spimquant tool command-line interface is a composition of the core 
(app-based) arguments and options, and the options related to Snakemake. 

As an end user, command line interface (CLI) is recommended as a way to specify 
the input/output data paths and configuration settings. The arguments passed 
in this way will be used to overwrite some of the default options defined in 
the config file `spimquant/config/snakebids.yml` and adapt to different datasets. 

Underlying SPIMquant, we use `snakemake` to manage workflow and config file. 
`spimquant` is a program that writes to the config file, and supports 
all Snakemake arguments. 

Steps to run a typical job looks like the following:
```bash
spimquant bids_dir output_dir {participant} [snakemake_args]
# now since we run spimquant, contents in config/snakebids.yml has changed
# and we can run snakemake with the updated configuration
snakemake result_file_path
```
All three arguments above to `spimquant` are not optional but necessary. 
The above first command will do the following:
1. Log the configuration used to run snakemake in the location output_dir, 
   including the snakebids.yml file and logs.
2. Pass bids_dir over to overwrite the bids_dir defined in the default 
   snakebids.yml. 
3. Start snakemake using the snakemake_args provided

A more concrete example:
```bash
spimquant test_bids_dir test_out participant
snakemake results/result_image.nii -n
```
Here we use -n so you can run this command but snakemake will only output a 
plan without generating any file. In an actual run, you can use `--cores all` 
to specify to use all cores and `--use-apptainer` to download and use existing 
container image for SPIMquant that have already installed all dependencies, 
without having to set up the environment yourself.

The core BIDS App arguments and app-specific options are listed below. 


```{argparse}
:ref: spimquant.run.get_parser
:prog: spimquant
```

## Snakemake-based command-line interface

Both core and Snakemake options are to be provided on the same command-line,  that 
is, any Snakemake CLI option can be used when running the app.

The Snakemake-based command-line options are detailed below. 

```{argparse}
:module: snakebids.snakemake_compat
:func: get_argument_parser
:prog: snakemake
```

