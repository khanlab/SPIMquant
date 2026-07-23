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


rule gaussian_biasfield:
    """simple bias field correction with gaussian"""
    input:
        spim=inputs["spim"].path,
    output:
        corrected=temp(
            bids_oz_out(
                root=work,
                datatype="seg",
                stain="{stain}",
                level="{level}",
                desc="correctedgaussian",
                suffix="SPIM.{ext}",
                **inputs["spim"].wildcards,
            ),
            group_jobs=True,
        ),
        biasfield=temp(
            bids_oz_out(
                root=work,
                datatype="seg",
                stain="{stain}",
                level="{level}",
                desc="gaussian",
                suffix="biasfield.{ext}",
                **inputs["spim"].wildcards,
            ),
            group_jobs=True,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=256000,
        disk_mb=2097152,
        runtime=15,
    params:
        proc_level=config["correction_level"],
        zarrnii_kwargs=zarrnii_in_kwargs,
    script:
        "../scripts/gaussian_biasfield.py"


rule n4_pre_quant:
    """Apply N4 bias field correction to SPIM images.
    
    Uses ANTs N4BiasFieldCorrection to correct intensity inhomogeneities within
    the brain mask. Outputs both the corrected image and the estimated bias field.
    """
    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        corrected=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="prequantN4",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="prequantN4",
            suffix="biasfield.nii.gz",
            **inputs["spim"].wildcards,
        ),
    params:
        iters=config["correction_n4_iters"],
        spline_spacing=config["correction_n4_spline_spacing"],
        shrink_level=1,  #shrink_level 1 since we use the correction_level to shrink
    threads: 16
    resources:
        mem_mb=32000,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -i {input.nii}"
        " -c [{params.iters}x{params.iters}x{params.iters}x{params.iters},0.0]"
        " -b [{params.spline_spacing},3]"
        " -s {params.shrink_level}"
        " -o [{output.corrected},{output.biasfield}]"
        " -x {input.mask}"
        " -d 3 -v "

rule calc_n4_rescaling:
    """ calculate the linear intensity rescaling,
     scale and offset, that n4 uses to rescale 
     intensities back to the input min and max.
     Does this by calculating the input min and max (within 
     the masked region), then calculating the min and max in the image
     input image divided by the bias field (also only within the masked 
     region). Then use these min/max values to find the scale and offset 
     (e.g. x*scale + offset) to apply to the divided image, so that the min and max
     become the input min and max. """
     input:
        uncorr=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="prequantN4",
            suffix="biasfield.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        scale_offset_params=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="prequantN4",
            suffix="scaleoffset.txt",
            **inputs["spim"].wildcards,
        ),
    script:
        "../scripts/calc_n4_rescaling.py" #TODO: make this script

rule encode_mask_in_bias_field:
    """Use -1 to encode the mask in the bias field"""
    input:
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level="{level}",
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="prequantN4",
            suffix="biasfield.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="prequantN4masked",
            suffix="biasfield.nii.gz",
            **inputs["spim"].wildcards,
        ),
    resources:
        mem_mb=32000,
        runtime=15,
    conda:
        "../envs/c3d.yaml"
    shell:
        "c3d {input.mask} {input.biasfield} -multiply " 
        " {input.mask} -replace 1 0 0 -1 "
        " -add -o {output.biasfield}"



ruleorder: n4_pre_quant > n4


rule n4_pre_quant_tune:
    """Grid-search over N4 parameters for parameter tuning QC.

    Runs N4BiasFieldCorrection for a single (spline_spacing, iters) combination
    encoded in the ``{n4_spline}`` and ``{n4_iters}`` wildcards.  Outputs are
    temporary; they are consumed by ``qc_n4_tune_png`` and ``qc_n4_tune_report``
    to produce the final PNG and HTML QC artefacts.
    """
    input:
        nii=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=str(config["correction_level"]),
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        mask=bids(
            root=root,
            datatype="micr",
            stain=stain_for_reg,
            level=str(config["correction_level"]),
            desc="brain",
            suffix="mask.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        corrected=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level=str(config["correction_level"]),
                desc="n4tuneSpline{n4_spline}Iters{n4_iters}",
                suffix="corrected.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
        biasfield=temp(
            bids(
                root=root,
                datatype="micr",
                stain="{stain}",
                level=str(config["correction_level"]),
                desc="n4tuneSpline{n4_spline}Iters{n4_iters}",
                suffix="biasfield.nii.gz",
                **inputs["spim"].wildcards,
            )
        ),
    params:
        iters="{n4_iters}",
        spline_spacing="{n4_spline}",
        shrink_level=1,
    threads: 16
    resources:
        mem_mb=32000,
        runtime=15,
    conda:
        "../envs/ants.yaml"
    shell:
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        "N4BiasFieldCorrection -i {input.nii}"
        " -c [{params.iters}x{params.iters}x{params.iters}x{params.iters},0.0]"
        " -b [{params.spline_spacing},3]"
        " -s {params.shrink_level}"
        " -o [{output.corrected},{output.biasfield}]"
        " -x {input.mask}"
        " -d 3 -v "


rule n4_biasfield:
    """N4 bias field correction with antspyx"""
    input:
        spim=inputs["spim"].path,
        biasfield=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=config["correction_level"],
            desc="prequantN4masked",
            suffix="biasfield.nii.gz",
            **inputs["spim"].wildcards,
        ),
        scale_offset_params=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level="{level}",
            desc="prequantN4",
            suffix="scaleoffset.txt",
            **inputs["spim"].wildcards,
        ),
    output:
        corrected=temp(
            bids_oz_out(
                root=work,
                datatype="seg",
                stain="{stain}",
                level="{level}",
                desc="correctedn4",
                suffix="SPIM.{ext}",
                **inputs["spim"].wildcards,
            ),
            group_jobs=True,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=500000 if config["dask_scheduler"] == "distributed" else 250000,
        runtime=180,
    params:
        proc_level=config["correction_level"],
        zarrnii_kwargs=zarrnii_in_kwargs,
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
        mask=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="otsu+k{k,[0-9]+}i{i,[0-9]+}",
            suffix="mask.{ext}",
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
    script:
        "../scripts/multiotsu.py"


rule gmmthresh:
    """Perform GMM tail-thresholding for segmentation.

Fits a Gaussian Mixture Model with n components to the log1p-transformed
global intensity histogram.  Components are sorted by mean intensity
ascending.  The threshold is defined as:

    threshold = expm1(mean_n + k * sigma_n)

where mean_n and sigma_n are the mean and standard deviation of the
highest-intensity component (nth component).  Outputs a QC figure
showing the histogram with GMM fits and threshold lines.

Method string format: gmm+n{n}k{k}, where k may use 'p' as decimal
separator (e.g. gmm+n2k2p5 → k=2.5).
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
        mask=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="gmm+n{n,[0-9]+}k{k,[0-9]+(?:p[0-9]+)?}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
        thresholds_png=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="gmm+n{n,[0-9]+}k{k,[0-9]+(?:p[0-9]+)?}",
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
        gmm_n=lambda wildcards: int(wildcards.n),
        gmm_k=lambda wildcards: float(wildcards.k.replace("p", ".")),
    script:
        "../scripts/gmmthresh.py"


rule threshold:
    """Apply simple intensity threshold for segmentation.

Creates binary mask by thresholding the corrected image at a specified intensity value.
Simpler alternative to multi-Otsu for cases where the threshold is known a priori.

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
        mask=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="th{threshold,[0-9]+}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    threads: 128 if config["dask_scheduler"] == "distributed" else 32
    resources:
        mem_mb=500000 if config["dask_scheduler"] == "distributed" else 250000,
        runtime=180,
    params:
        threshold=lambda wildcards: int(wildcards.threshold),
    script:
        "../scripts/threshold.py"


rule clean_segmentation:
    """Clean segmentation mask by removing edge artifacts and small objects.

Performs connected component analysis to identify and exclude objects that
extend too close to the image boundaries (likely artifacts). Creates both
a cleaned mask and an exclusion mask showing what was removed.
"""
    input:
        mask=bids_oz_in(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="{desc}",
            suffix="mask.{ext}",
            **inputs["spim"].wildcards,
        ),
    output:
        exclude_mask=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="{desc}+cleaned",
            suffix="excludemask.{ext}",
            **inputs["spim"].wildcards,
        ),
        cleaned_mask=bids_oz_out(
            root=root,
            datatype="seg",
            stain="{stain}",
            level="{level}",
            desc="{desc}+cleaned",
            suffix="mask.{ext}",
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
    script:
        "../scripts/clean_segmentation.py"
