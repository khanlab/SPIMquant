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
        Scheduler type: ``"threads"`` for the built-in threaded scheduler, or
        ``"distributed"`` for a ``dask.distributed.LocalCluster``.
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
        The distributed ``Client`` when *scheduler* is ``"distributed"``,
        otherwise ``None``.
    """
    if scheduler == "distributed":
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
            client.close()
            cluster.close()
    else:
        import dask

        dask.config.set(scheduler="threads", num_workers=threads)
        yield None
