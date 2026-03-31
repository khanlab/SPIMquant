"""Visual QC outputs for SPIMquant, generated at subject/channel level.

This module produces PNG quality-control figures covering:

1. Raw intensity / histogram plots
   - Per-channel linear and log-scale intensity histograms
   - Cumulative distribution and saturation/clip fraction

2. Segmentation overview figures
   - Slice montages (axial, coronal, sagittal) with field-fraction overlay
   - Max-intensity projection (MIP) with overlay

3. Spatial QC / coverage
   - Z-profile of mean signal intensity and segmented field fraction

4. Object-level summaries
   - Volume distribution, log-volume distribution, equivalent-radius
     distribution and summary statistics for detected objects

5. Per-ROI summaries (subject level)
   - Top-regions bar plots for field fraction, count, and density

All outputs are written to the ``qc`` datatype directory for each subject.
"""


rule qc_intensity_histogram:
    """Per-channel intensity histogram QC.

Reads the raw OME-Zarr at the registration downsampling level and
generates a four-panel figure: linear histogram, log-scale histogram,
cumulative distribution, and a summary-statistics panel including the
saturation/clip fraction (percentage of voxels at the maximum bin).
"""
    input:
        spim=inputs["spim"].path,
    output:
        png=bids(
            root=root,
            datatype="qc",
            stain="{stain}",
            suffix="histogram.png",
            **inputs["spim"].wildcards,
        ),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        level=config["registration_level"],
        hist_bins=500,
        hist_range=[0, 65535],
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/qc_intensity_histogram.py"


rule qc_segmentation_overview:
    """Segmentation overview slice montage QC.

Displays sample slices in axial, coronal, and sagittal orientations with
the segmentation field fraction overlaid on the SPIM background, and a
max-intensity projection column for each orientation.  Useful for quickly
identifying misregistration, segmentation artefacts, or coverage gaps.
"""
    input:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=config["registration_level"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        fieldfrac=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["registration_level"],
            desc="{desc}",
            suffix="fieldfrac.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        png=bids(
            root=root,
            datatype="qc",
            stain="{stain}",
            desc="{desc}",
            suffix="segslices.png",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=8000,
        runtime=15,
    script:
        "../scripts/qc_segmentation_overview.py"


rule qc_zprofile:
    """Z-profile QC: per-slice signal intensity and segmented fraction.

Plots the mean signal intensity (with ±1 SD band) and mean field fraction
across Z-slices.  Reveals depth-dependent artefacts such as striping,
illumination fall-off, or uneven staining.
"""
    input:
        spim=bids(
            root=root,
            datatype="micr",
            stain="{stain}",
            level=config["registration_level"],
            suffix="SPIM.nii.gz",
            **inputs["spim"].wildcards,
        ),
        fieldfrac=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["registration_level"],
            desc="{desc}",
            suffix="fieldfrac.nii.gz",
            **inputs["spim"].wildcards,
        ),
    output:
        png=bids(
            root=root,
            datatype="qc",
            stain="{stain}",
            desc="{desc}",
            suffix="zprofile.png",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10,
    script:
        "../scripts/qc_zprofile.py"


rule qc_objectstats:
    """Object-level statistics QC.

Loads the aggregated region-properties parquet (all stains combined) and
plots volume distribution, log-volume distribution, equivalent spherical
radius distribution, and a summary-statistics panel for the stain
specified by the ``{stain}`` wildcard.
"""
    input:
        regionprops=bids(
            root=root,
            datatype="tabular",
            desc="{desc}",
            space=config["template"],
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
    output:
        png=bids(
            root=root,
            datatype="qc",
            stain="{stain}",
            desc="{desc}",
            suffix="objectstats.png",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10,
    script:
        "../scripts/qc_objectstats.py"


rule qc_roi_summary:
    """Per-ROI summary QC: top-region bar plots for a single subject.

Reads the merged segmentation-statistics TSV (all stains) and the atlas
label table, then produces horizontal bar charts of the top brain regions
ranked by field fraction, object count, and density for every stain.
"""
    input:
        segstats=bids(
            root=root,
            datatype="tabular",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="mergedsegstats.tsv",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(
            root=root,
            template="{template}",
            seg="{seg}",
            suffix="dseg.tsv",
        ),
    output:
        png=bids(
            root=root,
            datatype="qc",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="roisummary.png",
            **inputs["spim"].wildcards,
        ),
    threads: 1
    resources:
        mem_mb=4000,
        runtime=10,
    script:
        "../scripts/qc_roi_summary.py"
