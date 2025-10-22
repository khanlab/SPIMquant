from dask.diagnostics import ProgressBar, Profiler, ResourceProfiler, CacheProfiler, visualize
from pathlib import Path
import cloudpickle

import dask
from dask.base import visualize

if __name__ == "__main__":

#    from dask.distributed import Client
#    client = Client()
    # or
 #   client = Client(processes=True)

    #from dask.distributed import LocalCluster
    #cluster = LocalCluster(processes=False, n_workers=1, threads_per_worker=128) # threads_per_worker=snakemake.threads)          # Fully-featured local Dask cluster
    #client = cluster.get_client()



    #from dask.distributed import performance_report




    dask.config.set(scheduler='threads')  # overwrite default with threaded scheduler

    def print_dask_diagnostics(graph, result=None, name="Dask task graph"):
        """Print diagnostics for a Dask graph and execution context."""
        if graph is not None:
            # Graph size
            num_tasks = len(graph.dask)
            print(f"\n=== {name} ===")
            print(f"Number of tasks in graph: {num_tasks}")
            print(f"Graph keys (first 5): {list(graph.dask.keys())[:5]}")
        else:
            print("\n(No Dask graph available)")

        # Scheduler info
        scheduler = dask.config.get("scheduler", "threads")
        num_workers = dask.config.get("num_workers", None)
        print(f"Using scheduler: {scheduler}")
        if num_workers:
            print(f"Configured workers: {num_workers}")

        # Thread count info
        import threading
        print(f"Active threads: {threading.active_count()}")

        # Example of computing the graph visualization (optional)
        # visualize(graph, filename="graph.svg")

        if result is not None:
            try:
                print(f"Result type: {type(result)}")
                if hasattr(result, 'shape'):
                    print(f"Result shape: {result.shape}")
            except Exception:
                pass
        print("==============================\n")


    from zarrnii.plugins import N4BiasFieldCorrection
    from zarrnii import ZarrNii

    hires_level = int(snakemake.wildcards.level)
    ds_level = int(snakemake.wildcards.dslevel)

    znimg = ZarrNii.from_ome_zarr(
        snakemake.input.spim,
        channel_labels=[snakemake.wildcards.stain],
        level=hires_level,
#        downsample_near_isotropic=True,
    )
    print(hires_level)
    print(ds_level)

 

    print(znimg)
    print(znimg.data)
    print_dask_diagnostics(znimg.data, name="Before bias field correction")

    #with performance_report(filename="dask-report.html"):
    with Profiler() as prof, ResourceProfiler() as rprof, CacheProfiler() as cprof:

        #with ProgressBar():

        print('perform bias field correction on downsampled, with upsampling')
        # apply bias field correction
        znimg_corrected = znimg.apply_scaled_processing(
            N4BiasFieldCorrection(sigma=5.0),
            downsample_factor=2**ds_level,
            chunksize=(1,20,20,20),
        )

        print_dask_diagnostics(znimg_corrected.data, name="before to_ome_zarr")

        # write to ome_zarr
        print('writing to ome zarr')
        znimg_corrected.to_ome_zarr(snakemake.output.corrected, max_layer=5)



    profiling_dir = Path(snakemake.output.profiling_dir)
    profiling_dir.mkdir(parents=True, exist_ok=True)



    # Save binary snapshot
    with open(profiling_dir / "profile.pkl", "wb") as f:
        cloudpickle.dump((prof, rprof, cprof), f)



