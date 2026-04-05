"""
Segmentation and quantification workflow for SPIMquant.

This module performs intensity correction, segmentation, and quantitative analysis of
pathology markers from SPIM microscopy data. It handles multi-channel data with different
stains (e.g., beta-amyloid, alpha-synuclein, Iba1) and produces per-region statistics.

Key workflow stages:
1. Intensity correction (Gaussian or N4 bias field correction)
2. Segmentation (threshold-based or multi-Otsu clustering)
3. Segmentation cleaning (remove edge artifacts and small objects)
4. Region properties computation (size, intensity, location of detected objects)
5. Coordinate transformation to template space
6. Colocalization analysis across channels
7. Atlas-based quantification (counts, density, field fraction per region)
8. Group-level statistical maps

Outputs include:
- Binary segmentation masks (OME-Zarr format)
- Region properties tables (Parquet format)
- Atlas-based statistics (TSV format)
- Volumetric density maps (NIfTI format)
"""

if config.get("direct_biasfield", False):

    rule gaussian_biasfield:
        """Gaussian bias field correction computed directly at the segmentation level (no downsample/upsample)"""
        input:
            spim=inputs["spim"].path,
        output:
            corrected=temp(
                directory(
                    bids(
                        root=work,
                        datatype="seg",
                        stain="{stain}",
                        level="{level}",
                        desc="correctedgaussian",
                        suffix="SPIM.ome.zarr",
                        **inputs["spim"].wildcards,
                    )
                ),
                group_jobs=True,
            ),
            biasfield=temp(
                directory(
                    bids(
                        root=work,
                        datatype="seg",
                        stain="{stain}",
                        level="{level}",
                        desc="gaussian",
                        suffix="biasfield.ome.zarr",
                        **inputs["spim"].wildcards,
                    )
                ),
                group_jobs=True,
            ),
        threads: 128 if config["dask_scheduler"] == "distributed" else 32
        resources:
            mem_mb=256000,
            disk_mb=2097152,
            runtime=15,
        params:
            zarrnii_kwargs={"orientation": config["orientation"]},
        script:
            "../scripts/gaussian_biasfield_direct.py"

    rule n4_biasfield:
        """N4 bias field correction computed directly at the segmentation level (no downsample/upsample)"""
        input:
            spim=inputs["spim"].path,
        output:
            corrected=temp(
                directory(
                    bids(
                        root=work,
                        datatype="seg",
                        stain="{stain}",
                        level="{level}",
                        desc="correctedn4",
                        suffix="SPIM.ome.zarr",
                        **inputs["spim"].wildcards,
                    )
                ),
                group_jobs=True,
            ),
        threads: 128 if config["dask_scheduler"] == "distributed" else 32
        resources:
            mem_mb=500000 if config["dask_scheduler"] == "distributed" else 250000,
            runtime=180,
        params:
            zarrnii_kwargs={"orientation": config["orientation"]},
            shrink_factor=16 if config["sloppy"] else 1,
            target_chunk_size=512,  #this sets the chunk size for this and downstream masks
        script:
            "../scripts/n4_biasfield_direct.py"

else:

    rule gaussian_biasfield:
        """simple bias field correction with gaussian"""
        input:
            spim=inputs["spim"].path,
        output:
            corrected=temp(
                directory(
                    bids(
                        root=work,
                        datatype="seg",
                        stain="{stain}",
                        level="{level}",
                        desc="correctedgaussian",
                        suffix="SPIM.ome.zarr",
                        **inputs["spim"].wildcards,
                    )
                ),
                group_jobs=True,
            ),
            biasfield=temp(
                directory(
                    bids(
                        root=work,
                        datatype="seg",
                        stain="{stain}",
                        level="{level}",
                        desc="gaussian",
                        suffix="biasfield.ome.zarr",
                        **inputs["spim"].wildcards,
                    )
                ),
                group_jobs=True,
            ),
        threads: 128 if config["dask_scheduler"] == "distributed" else 32
        resources:
            mem_mb=256000,
            disk_mb=2097152,
            runtime=15,
        params:
            proc_level=5,
            zarrnii_kwargs={"orientation": config["orientation"]},
        script:
            "../scripts/gaussian_biasfield.py"

    rule n4_biasfield:
        """N4 bias field correction with antspyx"""
        input:
            spim=inputs["spim"].path,
        output:
            corrected=temp(
                directory(
                    bids(
                        root=work,
                        datatype="seg",
                        stain="{stain}",
                        level="{level}",
                        desc="correctedn4",
                        suffix="SPIM.ome.zarr",
                        **inputs["spim"].wildcards,
                    )
                ),
                group_jobs=True,
            ),
        threads: 128 if config["dask_scheduler"] == "distributed" else 32
        resources:
            mem_mb=500000 if config["dask_scheduler"] == "distributed" else 250000,
            runtime=180,
        params:
            proc_level=5,
            zarrnii_kwargs={"orientation": config["orientation"]},
            shrink_factor=16 if config["sloppy"] else 1,
            target_chunk_size=512,  #this sets the chunk size for this and downstream masks
        script:
            "../scripts/n4_biasfield.py"


rule multiotsu:
    """Perform multi-Otsu thresholding for segmentation.

    Applies multi-level Otsu thresholding to identify intensity classes in the
    corrected image. The k parameter determines number of classes, and the i parameter
    selects which threshold index to use for creating the binary mask. Outputs a
    histogram visualization of the threshold selection.
    """
    input:
        corrected=bids(
            root=work,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="corrected{method}".format(method=config["correction_method"]),
            suffix="SPIM.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    output:
        mask=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="otsu+k{k,[0-9]+}i{i,[0-9]+}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
        thresholds_png=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="otsu+k{k,[0-9]+}i{i,[0-9]+}",
            suffix="thresholds.png",
            **inputs["spim"].wildcards,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=500000 if config["dask_scheduler"] == "distributed" else 250000,
        disk_mb=2097152,
        runtime=180,
    params:
        hist_bin_width=float(config["seg_hist_bin_width"]),
        hist_percentile_range=[float(x) for x in config["seg_hist_percentile_range"]],
        otsu_k=lambda wildcards: int(wildcards.k),
        otsu_threshold_index=lambda wildcards: int(wildcards.i),
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/multiotsu.py"


rule threshold:
    """Apply simple intensity threshold for segmentation.

    Creates binary mask by thresholding the corrected image at a specified intensity value.
    Simpler alternative to multi-Otsu for cases where the threshold is known a priori.
    """
    input:
        corrected=bids(
            root=work,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="corrected{method}".format(method=config["correction_method"]),
            suffix="SPIM.ome.zarr",
            **inputs["spim"].wildcards,
        ),
    output:
        mask=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="th{threshold,[0-9]+}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=500000 if config["dask_scheduler"] == "distributed" else 250000,
        runtime=180,
    params:
        threshold=lambda wildcards: int(wildcards.threshold),
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/threshold.py"


rule clean_segmentation:
    """Clean segmentation mask by removing edge artifacts and small objects.

    Performs connected component analysis to identify and exclude objects that
    extend too close to the image boundaries (likely artifacts). Creates both
    a cleaned mask and an exclusion mask showing what was removed.
    """
    input:
        mask=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    output:
        exclude_mask=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="{desc}+cleaned",
            suffix="excludemask.ozx",
            **inputs["spim"].wildcards,
        ),
        cleaned_mask=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="{desc}+cleaned",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=30,
    params:
        max_extent=0.15,
        proc_level=2,  #level at which to calculate conncomp
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/clean_segmentation.py"
