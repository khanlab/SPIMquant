rule run_vesselfm:
    """Run Vesselfm
    """
    input:
        "/nfs/trident3/lightsheet/prado/mouse_app_lecanemab_batch3/bids/sub-AS161F3/micr/sub-AS161F3_sample-brain_acq-imaris4x_SPIM.ome.zarr"
    params:
        
    output:
        "./new_res/"
    conda:
    "../envs/vesselfm.yaml"
    shell:
        "python -m vesselfm.cli --input-folder {input} --output-folder {output} --enable-dask-chunking --chunk-size 128 128 128 --downsample-level"