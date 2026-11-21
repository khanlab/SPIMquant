#!/usr/bin/env python3
"""Run a dask cluster with workers in Slurm jobs, for the jobs of one SPIMquant run.

Started by the workflow (``--dask_cluster slurm``) in a Slurm job of its own: it runs
the dask scheduler, submits the worker jobs with dask-jobqueue, and writes a scheduler
file that the workflow's jobs connect to (see ``dask_setup.get_dask_client``). It stops,
taking its worker jobs with it, when the stop file appears (written when the workflow
ends), when the workflow's heartbeat file has not been updated for
``--heartbeat-timeout`` seconds (the run was killed), or at the end of its time limit.
"""

import argparse
import os
import sys
import time
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scheduler-file", required=True)
    p.add_argument("--stop-file", required=True)
    p.add_argument("--log-dir", required=True)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--cores", type=int, required=True)
    p.add_argument("--mem-gb", type=int, required=True)
    p.add_argument("--processes", type=int, required=True)
    p.add_argument("--walltime", required=True)
    p.add_argument("--account", default=None)
    p.add_argument("--pythonpath", default="", help="prepended to workers' PYTHONPATH")
    p.add_argument("--heartbeat-file", required=True)
    p.add_argument("--heartbeat-timeout", type=float, default=900)
    p.add_argument(
        "--adaptive",
        action="store_true",
        help="scale between 0 and --jobs worker jobs with the load, rather than "
        "keeping --jobs",
    )
    p.add_argument(
        "--job-extra", action="append", default=[], help="extra #SBATCH directive"
    )
    return p.parse_args()


def _cluster_metrics(dask_scheduler):
    """Totals over the workers, from their latest heartbeats (runs on the scheduler)."""
    ws = list(dask_scheduler.workers.values())
    total = lambda key: sum(w.metrics.get(key, 0) or 0 for w in ws)
    return {
        "workers": len(ws),
        "threads": sum(w.nthreads for w in ws),
        "cpu": total("cpu"),
        "memory": total("memory"),
        "processing": sum(len(w.processing) for w in ws),
        "read": total("read_bytes"),
        "write": total("write_bytes"),
    }


def main():
    args = parse_args()

    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster

    # Workers read SPIM images through fsspec (e.g. mind://), and the workflow's
    # scripts send them functions from modules next to the scripts: workers need
    # the same environment (inherited from this job, as sbatch exports it) plus
    # the scripts directory on their path.
    prologue = [
        f"export PYTHONPATH={args.pythonpath}${{PYTHONPATH:+:$PYTHONPATH}}",
        # Return freed memory to the OS at once, as dask_setup does for local
        # clusters, so that unmanaged memory does not pause workers.
        "export MALLOC_TRIM_THRESHOLD_=0",
    ]

    cluster = SLURMCluster(
        cores=args.cores,
        processes=args.processes,
        memory=f"{args.mem_gb}GB",
        account=args.account,
        walltime=args.walltime,
        python=sys.executable,
        # spilled data goes to each node's local disk
        local_directory="$SLURM_TMPDIR",
        log_directory=args.log_dir,
        job_name="spimquant-dask-worker",
        job_script_prologue=prologue,
        job_extra_directives=args.job_extra,
        scheduler_options={"dashboard_address": ":8787"},
    )
    print(cluster.job_script(), flush=True)
    if args.adaptive:
        # keep idle workers for a while: the workflow's dask jobs come in bursts,
        # and a released node may take long to get back
        cluster.adapt(
            minimum_jobs=0, maximum_jobs=args.jobs, wait_count=60, interval="10s"
        )
    else:
        cluster.scale(jobs=args.jobs)

    client = Client(cluster)
    tmp = f"{args.scheduler_file}.tmp"
    client.write_scheduler_file(tmp)
    os.replace(tmp, args.scheduler_file)  # appears only once complete
    print(f"scheduler at {cluster.scheduler_address}", flush=True)
    print(f"dashboard at {cluster.dashboard_link} (host {os.uname().nodename})")
    print(f"requested {args.jobs} worker job(s)", flush=True)

    stop = Path(args.stop_file)
    heartbeat = Path(args.heartbeat_file)
    last_n = -1
    metrics = open(Path(args.log_dir) / "cluster_metrics.tsv", "a", buffering=1)
    metrics.write(
        "time\tworkers\tthreads\tcpu_pct_sum\tmemory_gb\tprocessing\tread_mbps\twrite_mbps\n"
    )
    try:
        while not stop.exists():
            age = time.time() - heartbeat.stat().st_mtime if heartbeat.exists() else 0
            if age > args.heartbeat_timeout:
                print(f"no heartbeat from the workflow for {age:.0f} s", flush=True)
                break
            if client.status != "running":
                print("lost the scheduler", flush=True)
                break
            n = len(client.nthreads())
            try:
                m = client.run_on_scheduler(_cluster_metrics)
                metrics.write(
                    f"{time.time():.0f}\t{m['workers']}\t{m['threads']}\t"
                    f"{m['cpu']:.0f}\t{m['memory'] / 1e9:.1f}\t{m['processing']}\t"
                    f"{m['read'] / 1e6:.1f}\t{m['write'] / 1e6:.1f}\n"
                )
            except Exception as e:  # metrics are a nicety
                print(f"metrics: {e}", flush=True)
            if n != last_n:
                print(f"{time.strftime('%H:%M:%S')} workers: {n}", flush=True)
                last_n = n
            time.sleep(10)
        else:
            print("stop file found", flush=True)
    finally:
        try:
            os.remove(args.scheduler_file)
        except OSError:
            pass
        client.close()
        cluster.close()  # cancels the worker jobs


if __name__ == "__main__":
    main()
