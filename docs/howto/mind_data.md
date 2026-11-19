# Processing Datasets on MIND

SPIMquant can process a BIDS dataset hosted on the [MIND](https://api.mind-western.ca) platform without
downloading it: give the dataset's `mind://` URL as the BIDS directory. Images are read from MIND as each job
needs them; outputs are written to a local (or cluster) directory.

!!! note "Early support"
    This uses development branches of snakebids (`remote-bids-dir`) and mind-client (`mind-fsspec`), installed
    by a separate pixi environment, `mind`. Once both are released, the default environment will do this too.

## Install

From a clone of SPIMquant:

```bash
BUILD_EDITABLE_PYTHON=false pixi install -e mind
```

`BUILD_EDITABLE_PYTHON=false` makes pixi install the snakebids branch as a regular package. Without it, pixi
installs it in editable mode, pointing into pixi's cache, so clearing that cache would break the environment
(reinstall with `BUILD_EDITABLE_PYTHON=false pixi reinstall -e mind snakebids` if that happens).

Run commands in this environment with `pixi run -e mind ...`, or enter it with `pixi shell -e mind`.

## Log in

```bash
pixi run -e mind mind-login
```

This prints an address: open it in any browser (it need not be on the machine you are logged in to), log in, and
paste back the whole address the browser lands on. The login is cached in `~/.config/fastapi_cfmm_client/`.

MIND's access tokens last only minutes, so for anything longer, keep them fresh in a second terminal (or `tmux`
window) for as long as SPIMquant runs:

```bash
pixi run -e mind mind-keep-token
```

## Run

Find the dataset's URL: `mind://` followed by the dataset's key prefix in the MIND store, e.g.
`mind://scratch/spimquant-test/test_mind_1x/bids`. (Until the MIND catalog lists dataset locations, ask whoever
uploaded the dataset.) Then run SPIMquant as usual, with the URL in place of the BIDS directory:

```bash
pixi run -e mind spimquant mind://scratch/spimquant-test/test_mind_1x/bids /path/to/output participant \
    --cores all
```

Any SPIMquant option works as with a local dataset. To check what will run first, add `-n` (a dry run, which
reads only the dataset's file listing and metadata).

The output's `dataset_description.json` records the source dataset's URL under `DatasetLinks`, so derivatives
stay traceable to the MIND dataset.

### On a Slurm cluster

Run SPIMquant from the login node with `--executor slurm` as usual. Jobs on compute nodes need the token too, so
keep a copy on a filesystem the nodes share, readable only by you, and point the jobs at it:

```bash
mkdir -m 700 -p ~/.mind
pixi run -e mind mind-keep-token --mirror ~/.mind/token.json     # in tmux, for the whole run

export MIND_TOKEN_FILE=~/.mind/token.json                         # in the shell that runs spimquant
pixi run -e mind spimquant mind://... /scratch/$USER/out participant --executor slurm -j 50
```

Slurm passes the shell's environment, including `MIND_TOKEN_FILE`, on to the jobs by default. Compute nodes need outbound HTTPS to the MIND server and to its
object store.

## How it works

- [mind-client](https://gitlab.com/cfmm/mind/mind_server/-/tree/main/client) registers `mind://` with
  fsspec and universal_pathlib, so pybids, snakebids, zarr and zarrnii read MIND objects by URL. Each read is
  signed by the MIND server, which redirects to the object store.
- snakebids indexes the remote dataset through pybids, and keeps `mind://` URLs as they are rather than treating
  them as local paths.
- Rules that read the raw SPIM image get its URL as a parameter rather than an input file, so Snakemake does not
  try to copy or track the remote image; local datasets still use it as a tracked input.

## Limitations

- **Read-only.** SPIMquant writes outputs locally; it does not upload them to MIND. Uploading derivatives back to
  MIND will be a separate step.
- **Throughput.** Each image chunk is a separate signed request, so whole-resolution steps (segmentation) are
  slower than on local storage. Run close to the data (e.g. on an Alliance cluster) for large datasets.
