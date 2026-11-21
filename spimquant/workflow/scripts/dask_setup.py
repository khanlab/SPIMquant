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
# Shared directory for temporary files, set by the workflow with --dask_cluster slurm
TMPDIR_ENV = "SPIMQUANT_DASK_TMPDIR"
# How long a job waits for the cluster: its scheduler and first worker are Slurm
# jobs that may still be queued when the job starts
WAIT_ENV = "SPIMQUANT_DASK_CLUSTER_WAIT_S"


@contextmanager
def connect_to_cluster():
    """Connect to the shared dask cluster of this run, waiting for it to start."""
    import os
    import time

    from dask.distributed import Client

    # Temporary files the workers write (zarr stores, including the one zarrnii
    # writes an .ozx through) must be on a filesystem they share with this job
    shared_tmp = os.environ.get(TMPDIR_ENV)
    if shared_tmp:
        import tempfile

        os.makedirs(shared_tmp, exist_ok=True)
        print(f"dask: temporary files in {shared_tmp} (was {tempfile.gettempdir()})")
        os.environ["TMPDIR"] = shared_tmp
        tempfile.tempdir = shared_tmp

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
        with _visible_stores(client):
            yield client
    finally:
        client.close()


@contextmanager
def _visible_stores(client):
    """Make zarr stores written by the cluster visible everywhere before they are read.

    On a shared filesystem such as Nibi's NFS-mounted /scratch, a chunk file that a
    worker on one node has written may for a short time not be visible on another
    node, even after the computation that wrote it has finished. Reading it there
    then finds no chunk, which zarr reads as the fill value (zeros): e.g. zarrnii,
    which builds each pyramid level from the level it stored before, made levels of
    zeros. Within one node (a LocalCluster) this does not happen.

    So, while the job is connected to the cluster:

    - zarr writes every chunk, including chunks of only the fill value, so that the
      chunk files a store must have are known from its shape;
    - every dask store into a local zarr array (``dask.array.store``, which
      ``to_zarr`` and zarrnii use) waits until all its chunk files can be opened on
      every node of the cluster and here, before returning.
    """
    import dask.array as da
    import dask.array.core as dac
    import zarr

    orig_store = dac.store

    def store(sources, targets, *args, **kwargs):
        result = orig_store(sources, targets, *args, **kwargs)
        if kwargs.get("compute", True):
            for target in targets if isinstance(targets, (list, tuple)) else [targets]:
                _wait_until_visible(client, target)
        return result

    with zarr.config.set({"array.write_empty_chunks": True}):
        dac.store = da.store = store
        try:
            yield
        finally:
            dac.store = da.store = orig_store


def _chunk_files(array):
    """The chunk (or shard) files of a zarr array in a local directory store, else None."""
    import itertools
    import os

    root = getattr(getattr(array, "store", None), "root", None)
    metadata = getattr(array, "metadata", None)
    if root is None or not hasattr(metadata, "encode_chunk_key"):
        return None
    prefix = os.path.join(str(root), array.path)
    return [
        os.path.join(prefix, metadata.encode_chunk_key(idx))
        for idx in itertools.product(*(range(n) for n in array.cdata_shape))
    ]


def _count_missing(paths):
    """How many of paths cannot be opened here, after refreshing their directories."""
    import os

    for d in {os.path.dirname(p) for p in paths}:
        try:
            os.listdir(d)  # revalidates the client's cached directory
        except OSError:
            pass
    missing = 0
    for p in paths:
        try:
            with open(p, "rb"):
                pass
        except OSError:
            missing += 1
    return missing


def _wait_until_visible(client, array, timeout=600):
    """Wait until every chunk file of array can be opened on each node and here."""
    import socket
    import time

    paths = _chunk_files(array)
    if not paths:
        return
    # one worker per node
    by_host = {}
    for addr, host in client.run(socket.gethostname).items():
        by_host.setdefault(host, addr)
    start = time.monotonic()
    while True:
        missing = client.run(_count_missing, paths, workers=list(by_host.values()))
        missing[socket.gethostname()] = _count_missing(paths)
        lagging = {h: n for h, n in missing.items() if n}
        if not lagging:
            waited = time.monotonic() - start
            if waited > 1:
                print(
                    f"dask: waited {waited:.0f} s for {len(paths)} chunk files to be visible"
                )
            return
        if time.monotonic() - start > timeout:
            raise TimeoutError(
                f"chunk files of {paths[0].rsplit('/c/', 1)[0]} still not visible after "
                f"{timeout} s: missing {lagging}"
            )
        time.sleep(1)


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
