"""Dask scheduler setup utilities for Snakemake workflow scripts."""

from contextlib import contextmanager


@contextmanager
def get_dask_client(
    scheduler, threads, threads_per_worker=2, mem_mb=None, reserve_mb=0
):
    """Context manager to set up the Dask scheduler.

    Parameters
    ----------
    scheduler : str
        Scheduler type: ``"threads"`` for the built-in threaded scheduler,
        ``"distributed"`` for a ``dask.distributed.LocalCluster``, or
        ``"cluster"`` to connect to the run's shared cluster
        (``--dask_cluster slurm``), in which case the other parameters do not
        apply.
    threads : int
        Total number of threads/cores available (from ``snakemake.threads``).
    threads_per_worker : int, optional
        Number of threads per worker when using the distributed scheduler.
        Ignored when *scheduler* is ``"threads"``.  Default is 2.
    mem_mb : int, optional
        Memory of the job in MB (``snakemake.resources.mem_mb``).  When given,
        the workers share ``mem_mb - reserve_mb`` evenly (capped at the memory
        dask detects, e.g. the SLURM cgroup limit).  Otherwise dask's
        ``"auto"`` limit is used, which gives each worker its share of the
        machine's memory *per core of the machine*: a job using half of a
        node's cores gets only half of its memory for the workers.
    reserve_mb : int, optional
        Memory in MB kept out of the workers' share, for the script itself
        (arrays it holds, results it collects) and the scheduler.

    Yields
    ------
    client : dask.distributed.Client or None
        The distributed ``Client`` when *scheduler* is ``"distributed"`` or
        ``"cluster"``, otherwise ``None``.
    """
    if scheduler == "cluster":
        # --dask_cluster slurm: connect to the run's shared cluster (started by
        # the workflow, scripts/dask_cluster_service.py) rather than starting one.
        # Workers are sized for the whole cluster, so threads_per_worker,
        # mem_mb and reserve_mb do not apply.
        with connect_to_cluster() as client, _performance_report():
            yield client
    elif scheduler == "distributed":
        import dask
        from dask.distributed import Client, LocalCluster

        # Return freed memory to the OS immediately. Without this, glibc
        # retains freed allocations as "unmanaged" worker memory, which
        # ratchets up until workers hit the 80% pause threshold and
        # deadlock (paused workers never drop back below the resume
        # threshold since the memory is never actually released).
        dask.config.set(
            {"distributed.nanny.pre-spawn-environ.MALLOC_TRIM_THRESHOLD_": 0}
        )

        n_workers = max(1, int(threads // threads_per_worker))
        memory_limit = "auto"
        if mem_mb:
            from distributed.system import MEMORY_LIMIT

            # snakemake's mem_mb goes to SLURM as --mem in MiB
            total = min(float(mem_mb) * 2**20, float(MEMORY_LIMIT))
            usable = max(total - float(reserve_mb) * 2**20, 0.25 * total)
            memory_limit = int(usable // n_workers)
        cluster = LocalCluster(
            n_workers=n_workers,
            threads_per_worker=threads_per_worker,
            memory_limit=memory_limit,
            dashboard_address=":8788",
        )
        client = Client(cluster)
        report = _performance_report()
        report.__enter__()
        print(cluster.dashboard_link)
        limits = [
            w["memory_limit"] for w in client.scheduler_info()["workers"].values()
        ]
        print(
            f"dask: {n_workers} workers x {threads_per_worker} threads, "
            f"{min(limits, default=0) / 2**30:.2f} GiB each",
            flush=True,
        )
        try:
            yield client
        finally:
            report.__exit__(None, None, None)
            client.close()
            cluster.close()
    else:
        import dask

        dask.config.set(scheduler="threads", num_workers=threads)
        yield None


# Set by the workflow with --dask_cluster slurm
SCHEDULER_FILE_ENV = "SPIMQUANT_DASK_SCHEDULER_FILE"
# How long a job waits for the cluster: its scheduler and first worker are Slurm
# jobs that may still be queued when the job starts
WAIT_ENV = "SPIMQUANT_DASK_CLUSTER_WAIT_S"


@contextmanager
def connect_to_cluster():
    """Connect to the shared dask cluster of this run, waiting for it to start."""
    import os
    import time

    from dask.distributed import Client

    scheduler_file = os.environ.get(SCHEDULER_FILE_ENV)
    if not scheduler_file:
        raise RuntimeError(
            f"dask scheduler 'cluster' needs {SCHEDULER_FILE_ENV} (set by "
            "--dask_cluster slurm)"
        )
    wait_s = float(os.environ.get(WAIT_ENV, 4 * 3600))
    start = time.monotonic()
    while not os.path.exists(scheduler_file):
        if time.monotonic() - start > wait_s:
            raise TimeoutError(
                f"no dask cluster after {wait_s:.0f} s: {scheduler_file} missing "
                "(is the spimquant-dask Slurm job queued or failed?)"
            )
        time.sleep(10)
    client = Client(scheduler_file=scheduler_file, timeout="120s")
    remaining = max(60, wait_s - (time.monotonic() - start))
    client.wait_for_workers(1, timeout=remaining)
    nthreads = client.nthreads()
    print(
        f"dask: connected to cluster {client.scheduler.address} "
        f"({client.dashboard_link}): {len(nthreads)} workers, "
        f"{sum(nthreads.values())} threads",
        flush=True,
    )
    try:
        yield client
    finally:
        client.close()


# Directory for dask performance reports (one HTML file per job), for benchmarking
REPORT_DIR_ENV = "SPIMQUANT_DASK_REPORT_DIR"


def _performance_report():
    """dask's performance report for the job, if SPIMQUANT_DASK_REPORT_DIR is set."""
    import contextlib
    import os
    import sys
    import time

    report_dir = os.environ.get(REPORT_DIR_ENV)
    if not report_dir:
        return contextlib.nullcontext()
    from dask.distributed import performance_report

    os.makedirs(report_dir, exist_ok=True)
    # snakemake runs a copy of the script named tmp<random>.<script>.py
    script = os.path.basename(sys.argv[0]).split(".")[-2] if sys.argv[0] else "dask"
    name = f"{script}_{time.strftime('%Y%m%d-%H%M%S')}_{os.getpid()}.html"
    return performance_report(filename=os.path.join(report_dir, name))
