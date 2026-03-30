# Workflow Rules Reference

This page lists every Snakemake rule in SPIMquant, grouped by the rule file that defines them.
Rules are the atomic units of computation in the workflow — each rule declares its inputs,
outputs, resources, and the shell command or Python script to run.

For a narrative description of what each *stage* does, see the
[Workflow Overview](workflow_overview.md).  For the files each stage produces, see
[Output Files](outputs.md).

---

## `import.smk` — Data Import

### `get_downsampled_nii`

Converts one OME-Zarr resolution level to a NIfTI file using
[zarrnii](https://github.com/khanlab/zarrnii).

- **Script**: `ome_zarr_to_nii.py`
- **Input**: SPIM OME-Zarr path (`inputs["spim"].path`)
- **Output**: `*_stain-{stain}_level-{level}_SPIM.nii.gz`
- **Threads**: 32 | **Memory**: 16 GB

### `import_template_anat`

Downloads or copies the template anatomy image from the URL or path in the config.

- **Output**: `tpl-{template}/tpl-{template}_anat.nii.gz`

### `import_template_spim`

Imports the SPIM-specific template if one is defined for the registration stain.

- **Output**: `tpl-{template}/tpl-{template}_stain-{stain}_anat.nii.gz`

### `import_mask`

Downloads or copies the template brain mask.

- **Output**: `tpl-{template}/tpl-{template}_desc-brain_mask.nii.gz`

### `import_dseg`

Downloads or copies the atlas discrete segmentation (parcellation).

- **Output**: `tpl-{template}/tpl-{template}_seg-{atlas}_dseg.nii.gz`

### `import_lut_tsv`

Downloads or copies the atlas label look-up table (TSV).

- **Output**: `tpl-{template}/tpl-{template}_seg-{atlas}_dseg.tsv`

### `copy_template_dseg_tsv`

Copies the atlas TSV to the output directory so downstream rules can find it.

- **Output**: `tpl-{template}/tpl-{template}_seg-{atlas}_dseg.tsv` (output copy)

### `generic_lut_bids_to_itksnap`

Converts a BIDS-format atlas TSV to ITK-SNAP colour table format.

- **Output**: `tpl-{template}/tpl-{template}_seg-{atlas}_dseg.txt`

---

## `masking.smk` — Brain Masking

### `pre_atropos`

Normalizes the SPIM intensity before Atropos GMM classification.

- **Input**: Downsampled SPIM NIfTI
- **Output**: Normalized SPIM NIfTI (temporary)

### `atropos_seg`

Runs ANTs Atropos GMM tissue classification with a spatial prior from the affine-registered
template mask.

- **Tool**: ANTs `Atropos`
- **Input**: Normalized SPIM + affine-warped template mask (prior)
- **Output**: Tissue probability maps + discrete segmentation (temporary)

### `post_atropos`

Post-processes the Atropos output (removes small objects, fills holes).

- **Output**: Cleaned tissue segmentation NIfTI (temporary)

### `init_affine_reg`

Quick affine registration of the downsampled SPIM to the template to provide a prior for
masking.

- **Tool**: `greedy`
- **Output**: Affine matrix (ITK text)

### `affine_transform_template_mask_to_subject`

Warps the template brain mask into subject space using the init affine transform.

- **Output**: `*_desc-brain_mask.nii.gz` (affine-transformed prior)

### `create_mask_from_gmm_and_prior`

Intersects the Atropos GMM result with the affine-warped template prior to produce the
final brain mask.

- **Output**: `*_stain-{stain}_level-{level}_desc-brain_mask.nii.gz`

---

## `templatereg.smk` — Template Registration

### `n4`

Applies ANTs N4 bias-field correction to the SPIM volume within the brain mask.

- **Tool**: `N4BiasFieldCorrection` (ANTs)
- **Output**: `*_desc-N4_SPIM.nii.gz`, `*_desc-N4_biasfield.nii.gz`
- **Memory**: 1.5 GB

### `apply_mask_to_corrected`

Multiplies the N4-corrected image by the brain mask.

- **Output**: `*_desc-N4masked_SPIM.nii.gz`

### `convert_ras_to_itk`

Converts the greedy transform (RAS convention) to ITK/ANTs convention for compatibility
with downstream ANTs tools.

- **Output**: ITK-format composite transform

### `affine_reg`

Performs rigid + affine registration of the masked, N4-corrected SPIM to the template.

- **Tool**: `greedy -a`
- **Output**: `*_desc-affine_xfm.txt`
- **Memory**: 16 GB | **Threads**: 8

### `init_affine_reg`

Generates an initial affine transform for the deformable registration.

- **Tool**: `greedy`
- **Output**: Initial affine matrix

### `deform_reg`

Performs diffeomorphic deformable registration of the affine-registered SPIM to the template.

- **Tool**: `greedy -d`
- **Output**: `*_desc-deform_xfm.nii.gz`, `*_desc-invdeform_xfm.nii.gz`
- **Memory**: ~64 GB | **Threads**: 16

### `compose_subject_to_template_warp`

Composes the affine and deformable transforms into a single composite warp field.

- **Tool**: `greedy` compose
- **Output**: `*_from-subject_to-{template}_desc-composite_xfm.nii.gz`

### `deform_spim_nii_to_template_nii`

Resamples the SPIM volume into template space using the composite warp.

- **Tool**: `antsApplyTransforms`
- **Output**: `*_space-{template}_SPIM.nii.gz`

### `deform_template_dseg_to_subject_nii`

Warps the atlas parcellation from template into subject space using the inverse deformation.

- **Tool**: `antsApplyTransforms` (nearest-neighbour interpolation)
- **Output**: `*_desc-deform_seg-{atlas}_dseg.nii.gz`

### `registration_qc_report`

Generates an HTML QC report overlaying the registered SPIM on the template.

- **Output**: `*_space-{template}_desc-reg_SPIM.html`
- **Memory**: 8 GB

---

## `segmentation.smk` — Pathology Segmentation

### `gaussian_biasfield`

Applies Gaussian smoothing-based bias-field correction to the raw OME-Zarr data.

- **Script**: `gaussian_biasfield.py`
- **Output**: Corrected OME-Zarr (`desc-correctedgaussian_SPIM.ome.zarr`, temporary)
- **Memory**: 256 GB | **Threads**: 128 (distributed) or 32 (thread)

### `n4_biasfield`

Applies ANTsPy N4 bias-field correction to the raw OME-Zarr data.

- **Script**: `n4_biasfield.py`
- **Output**: Corrected OME-Zarr (`desc-correctedn4_SPIM.ome.zarr`, temporary)
- **Memory**: 500 GB (distributed) or 250 GB (thread)

### `multiotsu`

Applies multi-Otsu threshold to segment pathology signal from the bias-corrected OME-Zarr.

- **Script**: `multiotsu.py`
- **Output**: Probability/threshold image OME-Zarr (temporary)

### `threshold`

Applies the Otsu-derived threshold to produce a binary candidate mask.

- **Script**: `threshold.py`
- **Output**: Binary mask OME-Zarr (temporary)

### `compute_filtered_regionprops`

Filters detected objects by minimum area and edge proximity; produces cleaned mask.

- **Script**: `compute_filtered_regionprops.py`
- **Output**: `*_desc-filtered_regionprops.parquet`, cleaned mask OME-Zarr

### `fieldfrac`

Downsamples the binary mask to registration resolution; computes field fraction (0–100%).

- **Script**: `fieldfrac.py`
- **Output**: `*_fieldfrac.nii.gz`

---

## `regionprops.smk` — Region Properties

### `map_regionprops_to_atlas_rois`

Assigns each detected object (from `regionprops.parquet`) to the atlas ROI that contains
its centroid.

- **Script**: `map_regionprops_to_atlas_rois.py`
- **Output**: `*_regionprops_with_roi.parquet` (temporary)

### `map_img_to_roi_tsv`

Converts per-object ROI assignments to a per-region summary (object count, field fraction,
volume).

- **Script**: `map_img_to_roi_tsv.py`
- **Output**: `*_desc-{atlas}_roi_tsv.parquet` (temporary)

### `transform_regionprops_to_template`

Applies the composite warp to centroid coordinates, mapping them to template space.

- **Script**: `transform_regionprops_to_template.py`
- **Output**: `*_space-{template}_regionprops.parquet`

### `aggregate_regionprops_across_stains`

Merges per-stain `regionprops.parquet` files and computes colocalization statistics.

- **Script**: `aggregate_regionprops_across_stains.py`
- **Output**: Multi-stain aggregated parquet

---

## `counts.smk` — Count Maps

### `counts_per_voxel`

Bins centroid coordinates into voxels at the registration resolution (subject space).

- **Script**: `counts_per_voxel.py`
- **Output**: `*_level-{level}_desc-filtered_count.nii.gz`

### `counts_per_voxel_template`

Bins template-space centroid coordinates into voxels (template space).

- **Script**: `counts_per_voxel_template.py`
- **Output**: `*_space-{template}_level-{level}_count.nii.gz`

---

## `fieldfrac.smk` — Field Fraction

### `deform_fieldfrac_nii_to_template_nii`

Warps the per-subject field fraction NIfTI into template space.

- **Tool**: `antsApplyTransforms` (Linear interpolation)
- **Output**: `*_space-{template}_fieldfrac.nii.gz`
- **Memory**: 32 GB | **Threads**: 32

---

## `segstats.smk` — Atlas Statistics

### `merge_into_segstats_tsv`

Merges the per-object ROI TSV and field fraction TSV into a single per-region statistics TSV.

- **Script**: `merge_into_segstats_tsv.py`
- **Output**: `*_stain-{stain}_desc-{atlas}_segstats.tsv` (temporary)

### `merge_into_colocsegstats_tsv`

Merges colocalization counts and field fraction into a colocalization statistics TSV.

- **Script**: `merge_into_segstats_tsv.py`
- **Output**: `*_desc-{atlas}_colocsegstats.tsv` (temporary)

### `merge_indiv_and_coloc_segstats_tsv`

Combines individual per-stain segstats with colocalization stats into a single wide-format TSV.

- **Script**: `merge_indiv_and_coloc_segstats_tsv.py`
- **Output**: `*_desc-{atlas}_mergedsegstats.tsv`

---

## `heatmaps.smk` — Heatmap Visualization

### `map_segstats_tsv_dseg_to_template_nii`

Paints each atlas region in template space with its segstats metric value.

- **Script**: `map_tsv_dseg_to_nii.py`
- **Output**: `*_space-{template}_desc-{atlas}_{metric}.nii.gz`
- **Memory**: 16 GB

### `map_segstats_tsv_dseg_to_subject_nii`

Paints each atlas region in subject space with its segstats metric value.

- **Script**: `map_tsv_dseg_to_nii.py`
- **Output**: `*_level-{level}_from-{template}_desc-{atlas}_{metric}.nii.gz`

### `deform_fieldfrac_nii_to_template_nii`

(Also in `fieldfrac.smk`) Warps field fraction to template space.

- **Tool**: `antsApplyTransforms`
- **Output**: `*_space-{template}_fieldfrac.nii.gz`

---

## `groupstats.smk` — Group Statistics

### `perform_group_stats`

Computes t-test and Cohen's *d* for each region × metric across contrast groups.

- **Script**: `perform_group_stats.py`
- **Input**: All subjects' `mergedsegstats.tsv` + `participants.tsv`
- **Output**: `group/*_groupstats.tsv`
- **Memory**: 1.5 GB

### `create_stats_heatmap`

Generates a PNG heatmap of the group statistical results.

- **Script**: `create_stats_heatmap.py`
- **Output**: `group/*_groupstats.png`
- **Memory**: 8 GB

### `map_groupstats_to_template_nii`

Paints the template atlas with group-level statistics (tstat/pval/cohensd).

- **Script**: `map_tsv_dseg_to_nii.py`
- **Output**: `group/*_{metric}_{stat}.nii.gz`
- **Memory**: 16 GB

### `concat_subj_parquet`

Concatenates per-subject `regionprops.parquet` files into a single group-level parquet.

- **Script**: `concat_subj_parquet.py`
- **Output**: `group/*_regionprops.parquet`

### `group_counts_per_voxel`

Creates a group-level voxel-wise density map in template space.

- **Script**: `counts_per_voxel_template.py`
- **Output**: `group/*_{stain}+count.nii.gz`
- **Memory**: 200 GB | **Threads**: 16

### `concat_subj_segstats_contrast` / `map_groupavg_segstats_to_template_nii`

Computes and maps per-contrast group-averaged regional metrics.

- **Scripts**: `concat_subj_segstats_contrast.py`, `map_tsv_dseg_to_nii.py`
- **Output**: `group/*_groupavgsegstats.tsv`, `group/*_groupavg.nii.gz`

---

## `vessels.smk` — Vessel Segmentation (optional)

### `import_vesselfm_model`

Copies the pre-trained VesselFM model weights to the local resource directory.

- **Output**: `resources/models/vesselfm.pt`

### `run_vesselfm`

Runs the VesselFM deep learning model to segment cerebrovascular structure.

- **Script**: `vesselfm.py`
- **Output**: `*_desc-vesselfm_mask.ozx`
- **GPU**: 1 | **Memory**: 64 GB

### `signed_distance_transform`

Computes a signed distance transform from the vessel mask (negative inside, positive outside).

- **Script**: `signed_distance_transform.py`
- **Output**: `*_dist.ozx`
- **Memory**: 256 GB | **Threads**: 128

---

## `patches.smk` — Patch Extraction (optional)

### `create_spim_patches`

Extracts 3D raw SPIM patches centred on atlas-region-sampled locations.

- **Script**: `create_patches.py`
- **Output**: `*_desc-raw_SPIM.patches/` directory
- **Memory**: 32 GB | **Threads**: 32

### `create_mask_patches`

Extracts 3D segmentation mask patches at the same locations as the SPIM patches.

- **Script**: `create_patches.py`
- **Output**: `*_desc-cleaned_mask.patches/` directory

### `create_corrected_spim_patches`

Extracts bias-corrected SPIM patches.

- **Script**: `create_patches.py`
- **Output**: `*_desc-corrected*_SPIM.patches/` directory

### `create_imaris_crops`

Extracts full-resolution bounding-box crops of atlas regions in Imaris format.

- **Script**: `create_imaris_crops.py`
- **Output**: `*_desc-crop_SPIM.imaris/` directory
- **Memory**: 32 GB | **Threads**: 32

---

## `preproc_mri.smk` — MRI Co-registration (optional)

This rule file is only included when `--register_to_mri` is set.  It co-registers an in-vivo
MRI (T2w or other structural scan) with the ex-vivo SPIM data.

Key rules:

| Rule | Description |
|------|-------------|
| `rigid_reg_mri_to_spim` | Rigid registration of MRI to SPIM |
| `deform_reg_mri_to_spim` | Deformable refinement of MRI–SPIM registration |
| `deform_spim_to_mri` | Warp SPIM into MRI space |
| `mri_registration_qc_report` | HTML QC report for MRI–SPIM registration |

See the [MRI Registration How-To](../howto/mri_registration.md) for details.

---

## Next Steps

- [Output Files Reference](outputs.md): Full list of output file paths and formats
- [Workflow Overview](workflow_overview.md): Narrative description of each stage
- [Configuration Options](config.md): How to tune rule behaviour