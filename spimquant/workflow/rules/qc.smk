"""Visual QC outputs for SPIMquant, generated at subject/channel level.

This module produces PNG quality-control figures covering:

1. Raw intensity / histogram plots
   - Per-channel linear and log-scale intensity histograms
   - Cumulative distribution and saturation/clip fraction

2. Segmentation overview figures
   - Slice montages (axial, coronal, sagittal) with field-fraction overlay,
     aspect ratio corrected from NIfTI voxel dimensions
   - Max-intensity projection (MIP) with overlay
   - Zoomed ROI montage: per-atlas-region crops with overlay detail

3. Spatial QC / coverage
   - Z-profile of mean signal intensity and segmented field fraction

4. Object-level summaries
   - Volume distribution, log-volume distribution, equivalent-radius
     distribution and summary statistics for detected objects

5. Per-ROI summaries (subject level)
   - Top-regions bar plots for field fraction, count, and density

6. Instance-centred animated GIFs
   - Per-stain: two GIFs (random order + sorted by radius) showing SPIM crops
     centred on each detected instance with all channels side-by-side, a
     hollow circle marker for the equivalent diameter, and atlas ROI annotation
   - Per-colocalization: same but centred on the midpoint of each colocalized
     object pair

Rules 2 and the ROI zoom are also generated for vessel segmentations.

All outputs are written to the ``qc`` datatype directory for each subject.
"""

# Whether to use the n4-corrected OME-Zarr as the background in zoom montage QC.
# Enabled via the --qc_roi_zoom_bg_n4 CLI option (only meaningful when
# correction_method=="n4").
_use_n4_bg = config.get("qc_roi_zoom_bg_n4", False) and (
    config.get("correction_method") == "n4"
)


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
    threads: 4
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

Loads SPIM data via ZarrNii (``downsample_near_isotropic=True``) and the
raw binary segmentation mask at the corresponding pyramid level.  Displays
sample slices in axial, coronal, and sagittal orientations with the mask
overlay, and a max-intensity projection column for each orientation.
Aspect ratio is corrected using voxel spacings from ``ZarrNii.get_zooms()``.
"""
    input:
        spim=inputs["spim"].path,
        mask=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ozx",
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
    threads: 4
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        level=config["registration_level"],
        mask_level=config["registration_level"] - config["segmentation_level"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/qc_segmentation_overview.py"


rule qc_vessels_overview:
    """Vessel segmentation overview slice montage QC.

Identical visualisation to ``qc_segmentation_overview`` but applied to the
vessel binary mask.  Loads data via ZarrNii with ``downsample_near_isotropic``
for isotropic display and physically correct aspect ratio.
"""
    input:
        spim=inputs["spim"].path,
        mask=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
    output:
        png=bids(
            root=root,
            datatype="qc",
            stain="{stain}",
            desc="{desc}",
            suffix="vesselslices.png",
            **inputs["spim"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        runtime=30,
    params:
        level=config["registration_level"],
        mask_level=config["registration_level"] - config["segmentation_level"],
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/qc_segmentation_overview.py"


rule qc_segmentation_roi_zoom:
    """Zoomed ROI montage QC for segmentation.

Crops the SPIM image and segmentation field-fraction mask to each atlas
region's bounding box (in subject space) and displays the best axial slice
with the field-fraction overlay.  Aspect ratio is corrected from NIfTI
voxel dimensions.  Provides detail-level visualisation of segmentation
quality within individual brain regions.

Two PNGs are produced: one with the mask overlay (``roimontage.png``) and
one without (``desc-{desc}nomask_roimontage.png``).
"""
    input:
        **(
            {
                "spim_n4": bids(
                    root=work,
                    datatype="seg",
                    stain="{stain}",
                    level=str(config["segmentation_level"]),
                    desc="correctedn4",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            }
            if _use_n4_bg
            else {}
        ),
        spim=inputs["spim"].path,
        mask=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
        dseg_nii=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
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
            stain="{stain}",
            desc="{desc}",
            suffix="roimontage.png",
            **inputs["spim"].wildcards,
        ),
        png_nomask=bids(
            root=root,
            datatype="qc",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}nomask",
            suffix="roimontage.png",
            **inputs["spim"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=15,
    params:
        max_rois=lambda wildcards: 25 if wildcards.seg == "coarse" else 100,
        n_cols=lambda wildcards: 5 if wildcards.seg == "coarse" else 10,
        patch_size=lambda wildcards: 2000 if wildcards.seg == "coarse" else 500,
        level=config["segmentation_level"],
        use_n4_bg=_use_n4_bg,
    script:
        "../scripts/qc_segmentation_roi_zoom.py"


rule qc_vessels_roi_zoom:
    """Zoomed ROI montage QC for vessel segmentation.

Identical to ``qc_segmentation_roi_zoom`` but applied to the vessel
binary mask.  Uses ZarrNii to load full-resolution data and
ZarrNiiAtlas for atlas-based ROI cropping.

Two PNGs are produced: one with the mask overlay (``vesselroimontage.png``)
and one without (``desc-{desc}nomask_vesselroimontage.png``).
"""
    input:
        **(
            {
                "spim_n4": bids(
                    root=work,
                    datatype="seg",
                    stain="{stain}",
                    level=str(config["segmentation_level"]),
                    desc="correctedn4",
                    suffix="SPIM.ome.zarr",
                    **inputs["spim"].wildcards,
                )
            }
            if _use_n4_bg
            else {}
        ),
        spim=inputs["spim"].path,
        mask=bids(
            root=root,
            datatype="vessels",
            stain="{stain}",
            level=config["segmentation_level"],
            desc="{desc}",
            suffix="mask.ozx",
            **inputs["spim"].wildcards,
        ),
        dseg_nii=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
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
            stain="{stain}",
            desc="{desc}",
            suffix="vesselroimontage.png",
            **inputs["spim"].wildcards,
        ),
        png_nomask=bids(
            root=root,
            datatype="qc",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}nomask",
            suffix="vesselroimontage.png",
            **inputs["spim"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=15,
    params:
        max_rois=lambda wildcards: 25 if wildcards.seg == "coarse" else 100,
        n_cols=lambda wildcards: 5 if wildcards.seg == "coarse" else 10,
        patch_size=lambda wildcards: 2000 if wildcards.seg == "coarse" else 500,
        level=config["segmentation_level"],
        use_n4_bg=_use_n4_bg,
    script:
        "../scripts/qc_segmentation_roi_zoom.py"


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


rule qc_otsu_threshold_sweep:
    """Threshold sweep QC HTML report for multiotsu segmentation.

Sweeps over a range of threshold values (spanning the configurable
percentile range of the bias-field corrected image) and generates 2D
crops at evenly-spaced positions in the image, one figure per threshold
value.  The otsu histogram PNG produced by the ``multiotsu`` rule is
embedded at the top of the report.  The resulting HTML report can be
visually assessed to select the optimal threshold before running the full
segmentation pipeline.

Only applicable when ``seg_method`` uses the ``otsu+k{}i{}`` pattern.
"""
    input:
        corrected=bids(
            root=work,
            datatype="seg",
            stain="{stain}",
            level=str(config["segmentation_level"]),
            desc="corrected{method}".format(method=config["correction_method"]),
            suffix="SPIM.ome.zarr",
            **inputs["spim"].wildcards,
        ),
        thresholds_png=bids(
            root=root,
            datatype="seg",
            stain="{stain}",
            level=str(config["segmentation_level"]),
            desc="{desc}",
            suffix="thresholds.png",
            **inputs["spim"].wildcards,
        ),
    output:
        html=bids(
            root=root,
            datatype="qc",
            stain="{stain}",
            desc="{desc}",
            suffix="otsuthreshqc.html",
            **inputs["spim"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=30,
    params:
        n_thresholds=10,
        n_crops=5,
        patch_size=300,
        level=config["segmentation_level"],
        hist_percentile_range=[float(x) for x in config["seg_hist_percentile_range"]],
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/qc_otsu_threshold_sweep.py"

rule qc_batch_otsu_threshold_sweep_raw:
    """Batch threshold sweep QC HTML report for multiotsu segmentation, on raw uncorrected images.

Uses the batch-wide Multi-Otsu thresholds (computed from the aggregated
histogram across all subjects in the dataset) to sweep over a range of
threshold values and display one mid-volume 2D crop per subject at each
threshold.  Subjects are shown as columns so that the user can immediately
see whether a given threshold generalises across the entire acquisition batch.

The aggregated histogram PNG is embedded at the top of the report.

Only applicable when ``seg_method`` uses the ``otsu+k{}i{}`` pattern.
"""
    input:
        zarrs=inputs["spim"].expand(),
        thresholds_png=bids(
            root=root,
            datatype="group",
            stain="{stain}",
            level=str(config["segmentation_level"]),
            desc="{desc}",
            suffix="thresholds.png",
        ),
    output:
        html=bids(
            root=root,
            datatype="group",
            stain="{stain}",
            desc="{desc}",
            suffix="batchotsuthreshqcraw.html",
        ),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=60,
    params:
        n_thresholds=10,
        patch_size=300,
        level=config["segmentation_level"],
        hist_percentile_range=[float(x) for x in config["seg_hist_percentile_range"]],
        zarrnii_kwargs=lambda wildcards: {"orientation": config["orientation"],"channel_labels":[wildcards.stain]},
    script:
        "../scripts/qc_batch_otsu_threshold_sweep.py"



rule qc_batch_otsu_threshold_sweep:
    """Batch threshold sweep QC HTML report for multiotsu segmentation.

Uses the batch-wide Multi-Otsu thresholds (computed from the aggregated
histogram across all subjects in the dataset) to sweep over a range of
threshold values and display one mid-volume 2D crop per subject at each
threshold.  Subjects are shown as columns so that the user can immediately
see whether a given threshold generalises across the entire acquisition batch.

The aggregated histogram PNG is embedded at the top of the report.

Only applicable when ``seg_method`` uses the ``otsu+k{}i{}`` pattern.
"""
    input:
        zarrs=inputs["spim"].expand(
            bids(
                root=work,
                datatype="seg",
                stain="{stain}",
                level=str(config["segmentation_level"]),
                desc="corrected{method}".format(method=config["correction_method"]),
                suffix="SPIM.ome.zarr",
                **inputs["spim"].wildcards,
            ),
            allow_missing=True,
        ),
        thresholds_png=bids(
            root=root,
            datatype="group",
            stain="{stain}",
            level=str(config["segmentation_level"]),
            desc="{desc}",
            suffix="thresholds.png",
        ),
    output:
        html=bids(
            root=root,
            datatype="group",
            stain="{stain}",
            desc="{desc}",
            suffix="batchotsuthreshqc.html",
        ),
    threads: 4
    resources:
        mem_mb=32000,
        runtime=60,
    params:
        n_thresholds=10,
        patch_size=300,
        level=config["segmentation_level"],
        hist_percentile_range=[float(x) for x in config["seg_hist_percentile_range"]],
        zarrnii_kwargs={"orientation": config["orientation"]},
    script:
        "../scripts/qc_batch_otsu_threshold_sweep.py"


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


rule qc_stain_instance_gif:
    """Instance-centred animated GIF QC for a single stain.

For every detected instance of the given stain, crops the raw SPIM data
(all channels side-by-side) centred on the instance position, overlays a
hollow circle whose diameter represents the equivalent spherical diameter,
and annotates with position, atlas ROI, and radius.  Two GIFs are produced:
one with a randomised frame order and one sorted by ascending equivalent
radius.

Inputs are the aggregated (all-stain) regionprops parquet in template space
(which retains subject-space ``pos_x/y/z`` columns) together with the atlas
parcellation for atlas-label lookup.
"""
    input:
        spim=inputs["spim"].path,
        instance_parquet=bids(
            root=root,
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        ),
        dseg_nii=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(
            root=root,
            template="{template}",
            seg="{seg}",
            suffix="dseg.tsv",
        ),
    output:
        gif_random=bids(
            root=root,
            datatype="qc",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="instancerandom.gif",
            **inputs["spim"].wildcards,
        ),
        gif_sorted=bids(
            root=root,
            datatype="qc",
            seg="{seg}",
            from_="{template}",
            stain="{stain}",
            desc="{desc}",
            suffix="instancesorted.gif",
            **inputs["spim"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        instance_type="stain",
        channels=stains,
        level=config["segmentation_level"],
        patch_size=config.get("instance_gif_patch_size", 100),
        max_instances=config.get("instance_gif_max_instances", 200),
        seed=42,
    script:
        "../scripts/qc_instance_gif.py"


rule qc_coloc_instance_gif:
    """Instance-centred animated GIF QC for colocalized object pairs.

Identical to ``qc_stain_instance_gif`` but operates on the colocalization
parquet: each frame shows the SPIM crop centred on the midpoint between the
two colocalized objects (subject-space ``pos_coloc_x/y/z``), with a circle
marker drawn at the average radius of the pair.
"""
    input:
        spim=inputs["spim"].path,
        instance_parquet=bids(
            root=root,
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="coloc.parquet",
            **inputs["spim"].wildcards,
        ),
        dseg_nii=bids(
            root=root,
            datatype="parc",
            seg="{seg}",
            level=config["registration_level"],
            from_="{template}",
            suffix="dseg.nii.gz",
            **inputs["spim"].wildcards,
        ),
        label_tsv=bids(
            root=root,
            template="{template}",
            seg="{seg}",
            suffix="dseg.tsv",
        ),
    output:
        gif_random=bids(
            root=root,
            datatype="qc",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="colocinstancerandom.gif",
            **inputs["spim"].wildcards,
        ),
        gif_sorted=bids(
            root=root,
            datatype="qc",
            seg="{seg}",
            from_="{template}",
            desc="{desc}",
            suffix="colocinstancesorted.gif",
            **inputs["spim"].wildcards,
        ),
    threads: 4
    resources:
        mem_mb=16000,
        runtime=60,
    params:
        instance_type="coloc",
        channels=stains,
        level=config["segmentation_level"],
        patch_size=config.get("instance_gif_patch_size", 100),
        max_instances=config.get("instance_gif_max_instances", 200),
        seed=42,
    script:
        "../scripts/qc_instance_gif.py"
