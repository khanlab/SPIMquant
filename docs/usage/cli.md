## Core command-line interface

The spimquant tool command-line interface is a composition of the core 
(app-based) arguments and options, and the options related to Snakemake. 

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

