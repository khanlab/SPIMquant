# Multi-node Dask Clusters on Slurm

By default, every job that processes a full-resolution image with Dask's distributed scheduler starts its
own Dask cluster (a `LocalCluster`) inside the job, so it can use at most the cores of one node. With
`--dask_cluster slurm`, SPIMquant instead starts one Dask cluster for the whole run, with its workers in
Slurm jobs (using [dask-jobqueue](https://jobqueue.dask.org)), and those jobs connect to it. A job can then
use several nodes, and the jobs running at the same time share the workers.

!!! note "Early support"
    This is new and has been tested on the Alliance cluster Nibi. Benchmarks are in the pull request that
    added it.

## Running

Run SPIMquant as usual (from a login node with `--executor slurm`, or inside a Slurm job), adding:

```bash
pixi run spimquant /path/to/bids /scratch/$USER/out participant \
    --dask_cluster slurm --dask_cluster_jobs 2 --dask_cluster_account rrg-mylab_cpu \
    --dask_cluster_walltime 12:00:00 ...
```

When the workflow starts, it submits a small Slurm job (`spimquant-dask`) that runs the Dask scheduler and
submits `--dask_cluster_jobs` worker jobs (`spimquant-dask-worker`). Jobs that use Dask wait until the
scheduler and at least one worker are running, so a queued worker job delays them rather than failing
them. When the workflow ends, the cluster stops and cancels its worker jobs; if the workflow is killed, the
cluster stops within about 15 minutes.

| Option | Default | |
|---|---|---|
| `--dask_cluster` | `local` | `slurm` for a shared, multi-node cluster |
| `--dask_cluster_jobs` | 1 | number of worker jobs (nodes, with the defaults below) |
| `--dask_cluster_cores` | 192 | cores per worker job (a whole Nibi node) |
| `--dask_cluster_mem_gb` | 700 | memory per worker job |
| `--dask_cluster_processes` | 96 | Dask worker processes per worker job (2 threads each) |
| `--dask_cluster_account` | `$SLURM_ACCOUNT` | Slurm account for the cluster's jobs |
| `--dask_cluster_walltime` | 12:00:00 | time limit of the cluster's jobs: cover the whole run |
| `--dask_cluster_adaptive` | off | start worker jobs only when there is work, and release them when idle |

The output directory must be on a filesystem shared by all nodes (e.g. `/scratch`), since the workers
write there.

!!! tip "Time limits and queue waits"
    The cluster's jobs, including the scheduler's, run for at most `--dask_cluster_walltime`, so set it to
    cover the run. On Alliance clusters, whole-node jobs of 3 hours or less go to a partition that starts
    much sooner: on Nibi, a whole node for 3 hours was estimated to start within 30 minutes, while one for
    4 hours was estimated to start 22 hours later.

### What changes for jobs that use Dask

- They request only 2 cores (`dask_cluster_client_threads` in the config): they build the task graph and
  collect results, while the cluster's workers do the work. They request 32 GB (`dask_cluster_client_mem_mb`).
- Their temporary files go to `.dask-cluster/tmp` in the output directory rather than the node's local
  disk, because workers on other nodes write to them.
- Per-rule Dask settings (threads per worker, worker memory) do not apply; the cluster's workers are sized
  by the options above.

Jobs that use the threaded scheduler (`--dask_scheduler threads`, or rules that always use it) are not
affected.

### Monitoring

The scheduler's log, the worker jobs' logs and a utilization table (`logs/cluster_metrics.tsv`: workers,
CPU, memory and tasks every 10 s) are in `.dask-cluster/` in the output directory. The scheduler log gives
the dashboard address (port 8787 on the scheduler's node); reach it with an SSH tunnel through the login
node, e.g. `ssh -L 8787:<node>:8787 nibi.alliancecan.ca`.

To get a Dask performance report (task stream, worker profiles) for every Dask job, set
`SPIMQUANT_DASK_REPORT_DIR=/path/to/reports` before running.

## When it helps

A multi-node cluster helps when a single job's full-resolution processing is limited by one node's cores,
or when several full-resolution jobs (several stains or subjects) run at the same time. Low-resolution
steps (registration, QC) are unaffected.
