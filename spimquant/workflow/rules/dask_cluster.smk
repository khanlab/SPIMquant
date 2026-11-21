"""
--dask_cluster slurm: one dask cluster, with workers in Slurm jobs, shared by the jobs of a run.

By default (--dask_cluster local) every job that uses the distributed scheduler starts a
LocalCluster on its own node, so a job can use at most one node. With --dask_cluster slurm,
the workflow starts a dask cluster when it starts (scripts/dask_cluster_service.py, in a
Slurm job of its own, which submits the worker jobs with dask-jobqueue) and stops it when it
ends. Jobs that use the distributed scheduler then connect to it (dask_setup.get_dask_client
with scheduler 'cluster'), so they can use several nodes, and several jobs share the workers.

What changes for those jobs:
- they run with --dask_cluster_client_threads threads, since the work is done by the
  cluster's workers, not the job itself;
- their tmpdir is on the output directory's filesystem, since workers on other nodes write
  the temporary zarr stores the jobs create there (and zarrnii writes .ozx outputs through a
  temporary directory, which the job's TMPDIR sets).

The cluster stops when the run ends (a stop file), or when the run's heartbeat file stops
being updated (the run was killed), or at --dask_cluster_walltime.
"""

if config.get("dask_cluster", "local") == "slurm":
    import shlex as _shlex
    import subprocess as _subprocess
    import sys as _sys
    import threading as _threading
    import time as _time

    if config["dask_scheduler"] not in ("distributed", "cluster"):
        raise ValueError("--dask_cluster slurm needs --dask_scheduler distributed")
    # dask_setup.get_dask_client connects to the shared cluster for this value
    config["dask_scheduler"] = "cluster"

    _dask_cluster_dir = os.path.abspath(os.path.join(root, ".dask-cluster"))
    _dask_scheduler_file = os.path.join(_dask_cluster_dir, "scheduler.json")
    _dask_stop_file = os.path.join(_dask_cluster_dir, "stop")
    _dask_heartbeat_file = os.path.join(_dask_cluster_dir, "heartbeat")
    _dask_jobid_file = os.path.join(_dask_cluster_dir, "slurm_jobid")
    _dask_tmpdir = os.path.join(_dask_cluster_dir, "tmp")

    # read by dask_setup in the jobs' scripts (each job parses this file again,
    # also when it runs in a Slurm job of its own)
    os.environ["SPIMQUANT_DASK_SCHEDULER_FILE"] = _dask_scheduler_file

    def _uses_dask_scheduler(rule):
        """Whether rule's script passes the --dask_scheduler setting to dask_setup."""
        if rule.script is None:
            return False
        path = rule.basedir.join(rule.script).get_path_or_uri(secret_free=True)
        try:
            with open(path) as f:
                return 'config["dask_scheduler"]' in f.read()
        except OSError:
            return False

    _dask_cluster_rules = [r for r in workflow.rules if _uses_dask_scheduler(r)]
    for _rule in _dask_cluster_rules:
        _rule.resources["_cores"] = int(config.get("dask_cluster_client_threads", 2))
        _rule.resources["tmpdir"] = _dask_tmpdir

    def _slurm_free_env():
        """The environment, without the variables of a Slurm job SPIMquant runs in.

        sbatch called from inside a job would otherwise pass them on (e.g. its memory
        per CPU, which conflicts with --mem).
        """
        return {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}

    def _slurm_time(value):
        """HH:MM:SS for Slurm. YAML reads an unquoted 4:00:00 as seconds (14400)."""
        if isinstance(value, int):
            return f"{value // 3600}:{value % 3600 // 60:02d}:{value % 60:02d}"
        return str(value)

    def _start_dask_cluster():
        os.makedirs(_dask_tmpdir, exist_ok=True)
        for f in (_dask_scheduler_file, _dask_stop_file):
            if os.path.exists(f):
                os.remove(f)
        account = config.get("dask_cluster_account") or os.environ.get("SLURM_ACCOUNT")
        walltime = _slurm_time(config["dask_cluster_walltime"])
        service = [
            _sys.executable,
            os.path.join(workflow.basedir, "scripts", "dask_cluster_service.py"),
            "--scheduler-file",
            _dask_scheduler_file,
            "--stop-file",
            _dask_stop_file,
            "--heartbeat-file",
            _dask_heartbeat_file,
            "--log-dir",
            os.path.join(_dask_cluster_dir, "logs"),
            "--jobs",
            str(config["dask_cluster_jobs"]),
            "--cores",
            str(config["dask_cluster_cores"]),
            "--mem-gb",
            str(config["dask_cluster_mem_gb"]),
            "--processes",
            str(config["dask_cluster_processes"]),
            "--walltime",
            walltime,
            "--pythonpath",
            os.path.join(workflow.basedir, "scripts"),
        ]
        if account:
            service += ["--account", account]
        if config.get("dask_cluster_adaptive", False):
            service += ["--adaptive"]
        sbatch = [
            "sbatch",
            "--parsable",
            "--job-name=spimquant-dask",
            f"--time={walltime}",
            "--cpus-per-task=4",
            "--mem=32G",
            f"--output={os.path.join(_dask_cluster_dir, 'scheduler-%j.log')}",
        ]
        if account:
            sbatch.append(f"--account={account}")
        sbatch += ["--wrap", _shlex.join(service)]
        jobid = (
            _subprocess.run(
                sbatch,
                check=True,
                capture_output=True,
                text=True,
                env=_slurm_free_env(),
            )
            .stdout.strip()
            .split(";")[0]
        )
        with open(_dask_jobid_file, "w") as f:
            f.write(jobid)
        logger.info(
            f"dask cluster: started Slurm job {jobid} for the scheduler, which "
            f"submits {config['dask_cluster_jobs']} worker job(s) of "
            f"{config['dask_cluster_cores']} cores; jobs that use dask wait for it. "
            f"Logs in {_dask_cluster_dir}"
        )

        def _heartbeat():
            while True:
                with open(_dask_heartbeat_file, "w") as f:
                    f.write(str(_time.time()))
                _time.sleep(60)

        _threading.Thread(target=_heartbeat, daemon=True).start()

    def _stop_dask_cluster():
        with open(_dask_stop_file, "w") as f:
            f.write("")
        jobid = (
            open(_dask_jobid_file).read().strip()
            if os.path.exists(_dask_jobid_file)
            else None
        )
        # the scheduler job checks for the stop file every 10 s and cancels its
        # workers; cancel it if it is still queued or does not stop
        for _ in range(12):
            if (
                jobid is None
                or not _subprocess.run(
                    ["squeue", "-h", "-j", jobid], capture_output=True, text=True
                ).stdout.strip()
            ):
                break
            _time.sleep(5)
        else:
            _subprocess.run(["scancel", jobid], env=_slurm_free_env())
        logger.info("dask cluster: stopped")

    _onstart_hooks.append(_start_dask_cluster)
    _onend_hooks.append(_stop_dask_cluster)
