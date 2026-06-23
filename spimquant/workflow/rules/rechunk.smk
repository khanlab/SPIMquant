"""Optional rechunking of input SPIM data to local work directory.

When --use-work-zarr is set, the rechunk_spim_to_work rule copies the full
multi-channel, multi-scale input SPIM OME-Zarr to the work directory with
configurable chunking.  The resulting ome.zarr directory is then used as the
input for downstream segmentation (gaussian_biasfield, n4_biasfield) and
vessel (run_vesselfm) rules instead of reading from the original source.

This is useful when the input data has suboptimal chunking, is stored on
slow network or cloud storage, or is in a format that requires decompression
on every random read (e.g. ome.zarr.zip).  The work_dir should be set to a
fast local disk (e.g. NVMe SSD) with sufficient storage capacity.
"""


rule rechunk_spim_to_work:
    """Copy input SPIM OME-Zarr to work directory with rechunking.

    Copies all channels and scale levels from the input OME-Zarr to the work
    directory.  The output is always an ome.zarr directory (not ozx), providing
    fast random access for downstream segmentation and vessel processing.

    Note: .ims (Imaris) format inputs are not supported by this rule.
    """
    input:
        spim=inputs["spim"].path,
    output:
        zarr=directory(work_zarr_path),
    threads: 32
    resources:
        mem_mb=64000,
        # 2 TiB: rechunking copies the full multi-scale, multi-channel zarr to
        # the work directory, which can be large for high-resolution datasets.
        disk_mb=2097152,
        runtime=120,
    params:
        chunks=config["zarr_store_chunks"],
    script:
        "../scripts/rechunk_spim.py"
