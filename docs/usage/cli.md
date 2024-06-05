The spimquant tool command-line interface is a composition of the core (app-based) 
arguments and options, and the options related to Snakemake. Both are provided on 
the same command-line, that is, any snakemake CLI option can be used when running spimquant.

## Core command-line interface

```{argparse}
:ref: spimquant.run.get_parser
:prog: spimquant
```

## Snakemake-based command-line interface

```{argparse}
:module: snakebids.snakemake_compat
:func: get_argument_parser
:prog: snakemake
```

