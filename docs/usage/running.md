# Running SPIMquant

```
spimquant BIDS_DIR OUTPUT_DIR ANALYSIS_LEVEL [options] [snakemake options]
```

From a clone, prefix with `pixi run` (or `pixi run -e gpu`). SPIMquant passes
every option it does not define to Snakemake, so the usual Snakemake options
work: `-n` (dry run), `--cores`, `--jobs`, `--executor`, `--until`,
`--forcerun`, `--rerun-incomplete`, `--keep-going`. The analysis levels and
SPIMquant's own options are in the [command-line reference](../reference/cli.md).

## On one machine

```bash
spimquant bids/ out/ participant -n            # list the jobs; nothing runs
spimquant bids/ out/ participant --cores all
spimquant bids/ out/ group --cores all
```

Jobs request memory (`mem_mb`) scaled from the image size. Snakemake only
counts memory against a limit if you set one, so on a workstation give it
the machine's memory: `--resources mem_mb=120000`.

A run that stops can be restarted with the same command; finished jobs are
not repeated. After a crash, add `--rerun-incomplete`. Each run's command is
recorded in `OUTPUT_DIR/code/snakebids/runs/`; `spimquant --history OUT`
lists them and `spimquant --replay OUT` runs the latest again.

## On a SLURM cluster

Run SPIMquant on the login node with the Snakemake SLURM executor; it submits
each job and waits:

```bash
spimquant bids/ out/ participant --executor slurm --jobs 100 --retries 1 --latency-wait 60
```

- The many short steps of each subject are submitted together, one job per
  stage: brain masking (`subj_mask`), registration (`subj_reg`), the region
  tables after segmentation (`subj_quant`) and, with two or more segmented
  stains, colocalization (`subj_coloc`).
- With `--retries 1`, a failed job is submitted once more with twice the
  memory and time for the steps whose resources grow with the image, and for
  ANTs registration. This also covers jobs that fail because a file written on
  another node is not yet visible on a shared filesystem.
- `--max-mem-mb` (default 450000) caps the memory any job requests, so a
  retried job still fits on a node.
- `--latency-wait 60` gives a shared filesystem a minute to show a file
  written on a compute node before Snakemake declares it missing (the
  built-in profile waits 10 s, which is too short on busy NFS).
- Run the login-node process in `tmux` or `screen`: it must stay alive
  until the last job finishes.

### GPU steps

LANTERN (`--seg_method lantern`, the Abeta default), vesselFM and FireANTs
need the `gpu` environment. pixi only installs it on a machine with a CUDA
driver, which login nodes usually lack, so tell it to assume one:

```bash
export CONDA_OVERRIDE_CUDA=12.8
pixi run -e gpu spimquant bids/ out/ participant --executor slurm --jobs 100 \
  --latency-wait 60 --retries 1
```

The login node only schedules; the GPU jobs run on GPU nodes with their real
driver. In the `default` environment, `--gpu auto` replaces default GPU
methods with CPU ones (Abeta: `gmm+n3k1`) and prints a note; a GPU method you
ask for explicitly is an error. vesselFM has no CPU replacement and runs on
the CPU instead, which is slow at full resolution; skip vessels with
`--no_vessels` if there is no GPU.

## Inside a job

Heavy steps (bias correction, segmentation, object measurement) run as dask
computations inside one Snakemake job, using its threads. By default each
starts a local distributed cluster (`--dask_scheduler distributed`); use
`threads` for a simpler, lower-overhead scheduler.

Intermediate OME-Zarr stores are written as single zip files (`.ozx`,
`--zarr-store-ext`), which keeps the file count low on shared filesystems.
`--work_dir` puts temporary files on fast local disk.

## Resources

Measured on one brain scanned at 1.625 × 1.625 × 2.75 µm (Abeta, CD31,
YoPro; 339 GB of Imaris), with the default settings, on the SLURM cluster
described above (GPU jobs on one A100). LANTERN and vesselFM read the raw
image, so no full-resolution bias correction runs. The whole participant level
took 3 h 40 min from submission. Steps that took under a minute are left out;
the brain-masking steps ran on a 12-core workstation. See
[Example outputs](../examples/index.md) for this run's outputs.

<!-- RESOURCES-TABLE: from docs/gen/collect_examples.py on the example run; see AGENTS.md -->
| Rule | Jobs | Longest job | Peak memory |
|---|---|---|---|
| `run_vesselfm` | 1 | 2.9 h | 31.3 GB |
| `signed_distance_transform` | 1 | 31 min | 92.8 GB |
| `run_lantern_plaques` | 1 | 22 min | 35.5 GB |
| `compute_filtered_regionprops` | 1 | 10 min | 31.2 GB |
| `binarize_lantern_plaques` | 1 | 9 min | 48.3 GB |
| `sample_sdt_at_objects` | 1 | 6 min | 0.7 GB |
| `qc_vessels_roi_zoom` | 4 | 5 min | 3.4 GB |
| `qc_segmentation_roi_zoom` | 4 | 5 min | 3.4 GB |
| `deform_reg` | 1 | 4 min | 12.5 GB |
| `init_affine_reg` | 1 | 2 min | 3.1 GB |
| `compose_subject_to_template_warp` | 1 | 2 min | 5.8 GB |
| `affine_reg` | 1 | 1 min | 3.4 GB |
| `n4_pre_quant` | 1 | 1 min | 0.1 GB |
| `registration_qc_report` | 1 | 1 min | 4.5 GB |
| `atropos_seg` | 1 | 1 min | 0.2 GB |

Total over all jobs: 5.3 h of wall time (jobs run in parallel).
<!-- /RESOURCES-TABLE -->

To measure your own data, add `--benchmark`: each job's run time, peak memory
and CPU use are written to `OUTPUT_DIR/benchmarks/`, and

```bash
pixi run -e dev python tools/summarize_benchmarks.py OUT [OUT2 ...]
```

compares what each rule requested with what it used.

## Checking a run

- `spimquant ... -n` lists the jobs that would run and why.
- `spimquant ... --rulegraph | dot -Tsvg > rules.svg` draws the rules of your
  run; the [workflow graphs](../reference/workflow.md) show the defaults.
- Logs of each job are under `OUTPUT_DIR/logs/`.
- Before processing a new dataset, the [dataset QC report](../howto/dataset_qc.md)
  checks every scan's orientation, stains and coverage in minutes.
