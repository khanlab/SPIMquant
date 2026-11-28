# Contributing

Bug reports, questions and pull requests are welcome at
<https://github.com/khanlab/SPIMquant>.

## Setup

```bash
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant
pixi install -e dev
```

## Checks

CI runs these on every push; run them before opening a pull request:

```bash
pixi run -e dev quality_check        # black, isort, snakefmt (quality_fix to fix)
pixi run -e dev pytest tests/        # unit tests of the workflow scripts
pixi run spimquant tests/bids_ds /tmp/out participant -n   # the workflow builds
```

The dry runs in `.github/workflows/test.yml` cover the main options (MRI
registration, cohort templates, fibres, job groups for SLURM). A pull request
that adds an option should add a dry run that uses it.

## Code layout

| Path | Contents |
|---|---|
| `spimquant/run.py` | The BIDS App entry point (snakebids). |
| `spimquant/config/snakebids.yml` | CLI options (`parse_args`), defaults, templates. |
| `spimquant/workflow/Snakefile` | Configuration logic and the target rules of each analysis level. |
| `spimquant/workflow/rules/*.smk` | Rules, by stage. |
| `spimquant/workflow/scripts/` | Python scripts run by the rules. |
| `tests/` | Unit tests and the small test datasets. |
| `docs/` | This documentation. |

Rules write BIDS-style names with `bids()`; files kept for the user go under
the output directory (`root`), and temporary ones are marked `temp()`.
Resources (`mem_mb`, `runtime`) of heavy rules scale with the image size; see
[Running SPIMquant](usage/running.md#resources) for how to measure them.

## Documentation

Most reference pages are generated from the code when the docs are built:
the command-line reference from the parser, the templates and configuration
keys from `spimquant/config/snakebids.yml`, and the output files and
workflow graphs from dry runs of the test dataset. What you need to do when
you change the workflow:

| You changed | Do |
|---|---|
| A CLI option or its help text | Nothing: the [command-line reference](reference/cli.md) is built from it. Write help text for users, with units and the default. |
| A config key without a CLI option | Put a comment above it in `snakebids.yml`; it becomes the description in the [configuration reference](reference/config.md). |
| Rules or output files | `pixi run docs-introspect`, then commit `docs/generated/`. CI fails if it is out of date. |
| A new output suffix, folder, template or analysis level | Describe it in `docs/gen/descriptions.yml`. `pixi run -e docs docs-check` lists what is missing; CI runs it too. |
| A new optional output | Add a scenario that turns it on to `docs/gen/scenarios.yml`, then as above. |
| QC reports or figures | Refresh the [example outputs](examples/index.md) from a full-resolution run: see `AGENTS.md`. |
| Behaviour described in a guide | Update the guide (`docs/howto/`, `docs/usage/`). |

Build and preview the site with:

```bash
pixi run -e docs docs-serve    # http://127.0.0.1:8000
pixi run -e docs docs-build    # strict: fails on broken links
```

Write as in the existing pages: say what an option does, to which files, with
the numbers a user needs. No marketing adjectives, no restating headings.
