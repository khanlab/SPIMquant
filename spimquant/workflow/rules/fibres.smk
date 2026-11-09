"""
Fibre orientation analysis (--fibres), after Foa3D (Sorelli et al. 2023).

A multiscale Frangi (Hessian) filter of the bias-corrected image gives the
fibre axis at each voxel; fibre voxels are reduced to additive per-super-voxel
moments (counts, orientation tensor and spherical harmonic sums), from which
orientation distribution functions (ODFs), dispersion and anisotropy maps are
derived.  See scripts/fibre_orientation.py.
"""


rule fibre_moments:
    """Multiscale Frangi filter and per-super-voxel fibre orientation moments.

    The Frangi background parameter (gamma, per scale) and the fibre threshold
    are global, estimated from voxels sampled in tissue-rich chunks, and saved
    to the params JSON with a QC histogram.  The moments are written as an
    OME-Zarr with one channel per moment, on a grid of super-voxels of about
    fibre_block_um.
    """
    input:
        corrected=bids_oz_in(
            root=work,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="corrected{method}".format(method=config["correction_method"]),
            suffix="SPIM.{ext}",
            **inputs["spim"].wildcards,
        ),
    output:
        moments=bids_oz_out(
            root=root,
            datatype="fibres",
            stain="{stain}",
            level="{level}",
            desc="frangi",
            suffix="moments.{ext}",
            **inputs["spim"].wildcards,
        ),
        params_json=bids(
            root=root,
            datatype="fibres",
            stain="{stain}",
            level="{level}",
            desc="frangi",
            suffix="params.json",
            **inputs["spim"].wildcards,
        ),
        qc_png=bids(
            root=root,
            datatype="fibres",
            stain="{stain}",
            level="{level}",
            desc="frangi",
            suffix="thresholds.png",
            **inputs["spim"].wildcards,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        # dask tasks are bounded by fibre_chunk_size (~0.45 GB per 128^3 task);
        # the moments are 1/512 of the voxels or less, so little else grows
        # with the image.
        mem_mb=scaled_mem_mb(64000, 0),
        disk_mb=2097152,
        runtime=scaled_runtime(30, 0.25),
    params:
        dask_threads_per_worker=4,
        dask_reserve_mb=8000,
        scales_um=config["fibre_scales_um"],
        block_um=config["fibre_block_um"],
        lmax=config["fibre_lmax"],
        alpha=config["fibre_alpha"],
        beta=config["fibre_beta"],
        gamma_quantile=config["fibre_gamma_quantile"],
        gamma_factor=config["fibre_gamma_factor"],
        threshold_method=config["fibre_threshold_method"],
        stat_blocks=config["fibre_stat_blocks"],
        stat_samples_per_block=config["fibre_stat_samples_per_block"],
        stat_level_offset=config["fibre_stat_level_offset"],
        edge_min_weight=config["fibre_edge_min_weight"],
        chunk_size=config["fibre_chunk_size"],
    script:
        "../scripts/fibre_moments.py"


rule fibre_metrics:
    """Fibre orientation maps (subject space) from the moments.

    Scalar maps (fibre density, mean Frangi response, mean scale, dispersion
    indices odi/odipri/odisec/odianis, orientation tensor FA, ODF GFA), the
    principal fibre direction (dir, world x/y/z) and its colour map (dec), and
    the ODF as real SH coefficients in the MRtrix3 convention (odf).
    """
    input:
        moments=bids_oz_in(
            root=root,
            datatype="fibres",
            stain="{stain}",
            level="{level}",
            desc="frangi",
            suffix="moments.{ext}",
            **inputs["spim"].wildcards,
        ),
    output:
        **{
            metric: bids(
                root=root,
                datatype="fibres",
                stain="{stain}",
                level="{level}",
                desc="frangi",
                suffix=f"{metric}.nii.gz",
                **inputs["spim"].wildcards,
            )
            for metric in FIBRE_SCALAR_METRICS + ["dir", "dec", "odf"]
        },
    group:
        "subj_quant"
    threads: 4
    resources:
        mem_mb=scaled_mem_mb(4000, 2),
        runtime=10,
    script:
        "../scripts/fibre_metrics.py"
