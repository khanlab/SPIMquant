# Multi-node dask on Nibi (reference, not merged)

This branch (`dask-slurm-cluster`) adds two experimental options, kept for reference rather
than merged:

- `--dask_cluster slurm`: one dask cluster for the whole run, with workers in Slurm jobs
  (dask-jobqueue), shared by the jobs that use dask, possibly across nodes.
- `--zarr_store_root`: OME-Zarr outputs written to object storage (e.g. MIND), tracked by
  local `.zarr.url` pointer files.

`report.html` benchmarks both on Nibi (2026-11-21), with a dataset read from MIND: bias-field
correction, GMM thresholding and region properties of three stains at full resolution. The
`run_*.json` files hold the data behind it.

Findings:

- Both work and give the same results as the default (per-job LocalCluster, local `.ozx`).
- Neither made the run faster. The steps were limited by reading from MIND, not by cores
  (workers averaged 5–11 cores busy of 192–576), and the cluster added Slurm queue waits
  for its worker jobs and used 1.3–3.5× the core-hours.
- Writing outputs to MIND made later steps 2–5× slower, as they read them back from MIND.
- With workers on several nodes, chunk files written on one node were not always visible
  on another right afterwards (Nibi's NFS /scratch), which silently gave pyramid levels of
  zeros. The branch works around it with a barrier after each store.

When it might be worth revisiting: a single image too large for one node, in a step that
is limited by compute or memory rather than by reading its input.
