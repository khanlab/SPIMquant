# Command Line Interface Reference

This page provides auto-generated documentation for all SPIMquant command-line options.

## Overview

```bash
spimquant <bids_dir> <output_dir> <analysis_level> [options]
```

## Positional Arguments

### `bids_dir`

The directory with the input dataset formatted according to the BIDS standard.
 (Path)



### `output_dir`

The directory where the output files should be stored. If you are running group level analysis this folder should be prepopulated with the results of the participant level analysis.
 (Path)



### `analysis_level`

Level of the analysis that will be performed.
 (str)
 (choices: participant, group)


## Optional Arguments

### `--work_dir`, `--work-dir`

Local path to use for temporary files
 (str)


### `--template`

Template to use for SPIM registration
 (str)
 (choices: ABAv3, DSURQE, gubra, MBMv3, turone)
 (default: ABAv3)


### `--template-mri`, `--template_mri`

Template to use for MRI registration to obtain brain mask
 (str)
 (choices: MouseIn, DSURQE, MBMv3, turone)
 (default: MouseIn)


### `--atlas_segs`, `--atlas-segs`

Atlas segmentations to use with the chosen template (default: use them all)
 (str)
 (accepts one or more values)


### `--patch-atlas-segs`, `--patch_atlas_segs`

Atlas segmentations to use for extracting patches (default: roi22)
 (str)
 (default: roi22)
 (accepts one or more values)


### `--template-negative-mask`, `--template_negative_mask`

Negative mask, in the template space, to highlight regions to avoid
 (Path)
 (default: placeholder)


### `--template-crop`, `--template_crop`

Crop template along X-axis to retain specific hemisphere for registration (default: None)
 (str)
 (choices: left, right)


### `--stains-for-reg`, `--stains_for_reg`

Possible stains to use for registration (will choose first available, in order)  (default: ['PI', 'YOPRO', 'YoPro', 'AutoF', 'autof'])
 (str)
 (default: PI, YOPRO, YoPro, AutoF, autof)
 (accepts one or more values)


### `--stains_for_seg`, `--stains-for-seg`

List of stains to use for segmentation and quantification  (default: ['abeta', 'Abeta', 'BetaAmyloid', 'AlphaSynuclein', 'Iba1', 'ChAT'])
 (str)
 (default: abeta, Abeta, BetaAmyloid, AlphaSynuclein, Iba1, ChAT)
 (accepts one or more values)


### `--registration-level`, `--registration_level`

Downsampling level to use for registration (level 0 is full res, level 1 is 50% size, ...) (default: 5)
 (str)
 (default: 5)


### `--segmentation-level`, `--segmentation_level`

Downsampling level to use for segmentation (level 0 is full res, level 1 is 50% size, ...) (default: 0)
 (str)
 (default: 0)


### `--no-segmentation`, `--no_segmentation`

Skip segmentation and quantification of stains, i.e. perform registration only (default: False)
 (default: False)
 (flag)


### `--correction_method`, `--correction-method`

Method to use for intensity non-uniformity correction, prior to performing segmentation (default: n4)
 (str)
 (choices: gaussian, n4)
 (default: n4)


### `--seg_method`, `--seg-method`

Method to use for microscopy segmentation (e.g. plaques, protein deposits, cells) applied to 'stains_for_seg' channels, and used to calculate field fractions.
 (str)
 (default: otsu+k3i2, th900)
 (accepts one or more values)


### `--seg_hist_range`, `--seg-hist-range`

Range of intensities to use for histogram calculation in multiotsu segmentation. Only applicable when seg_method is otsu+k{}i{}. Specify 2 numbers, for min and max values. (default: [0, 1000])
 (str)
 (default: 0, 1000)
 (accepts 2 values)


### `--seg_hist_bins`, `--seg-hist-bins`

Number of bins to use for histogram calculation in multiotsu segmentation. Only applicable when seg_method is otsu+k{}i{}. (default: 1000)
 (str)
 (default: 1000)


### `--register_to_mri`, `--register-to-mri`

Register the lightsheet data directly to a corresponding MRI (from the BIDS dataset)
 (default: False)
 (flag)


### `--orientation`

Forcefully set the orientation of the input OME Zarr. Only use this if the orientation in the zarr metadata is incorrect and not easily changed. Note, this orientation string is defined in XYZ, with the letter of the increasing direction. Blaze data stitched with SPIMprep should be RPI, Imaris stitched data should be RAI, and Lifecanvas stitched data should be RAS.
 (str)


### `--sloppy`

Use low-quality parameters for speed (USE FOR TESTING ONLY)
 (default: False)
 (flag)


### `--skip_bids_validation`, `--skip-bids-validation`

Skip validation of BIDS dataset. BIDS validation is performed by default using the bids-validator plugin (if installed/enabled) or with the pybids validator implementation (if bids-validator is not installed/enabled).
 (default: False)
 (flag)


### `--patch-size`, `--patch_size`

Size of patches to extract in voxels (x, y, z) (default: [256, 256, 256])
 (int)
 (default: 256, 256, 256)
 (accepts 3 values)


### `--n-patches-per-label`, `--n_patches_per_label`

Number of patches to extract per atlas label (default: 5)
 (int)
 (default: 5)


### `--patch_labels`, `--patch-labels`

List of atlas label names, abbreviations, or indices to extract patches from. If not specified, patches are extracted from all labels.
 (str)
 (accepts one or more values)


### `--patch-seed`, `--patch_seed`

Random seed for reproducible patch sampling (default: 42)
 (int)
 (default: 42)


### `--crop_labels`, `--crop-labels`

List of atlas label names, abbreviations, or indices to extract as Imaris crops. If not specified, all labels are cropped.
 (str)
 (accepts one or more values)


### `--crop-atlas-segs`, `--crop_atlas_segs`

Atlas segmentations to use for extracting Imaris crops (default: roi22)
 (str)
 (default: roi22)
 (accepts one or more values)


### `--contrast-column`, `--contrast_column`

Column name in participants.tsv to use for defining group contrasts (e.g., 'treatment', 'genotype'). Required for group-level statistical analysis.
 (str)


### `--contrast-values`, `--contrast_values`

Two group values for contrast comparison (e.g., 'control' 'drug'). Used with --contrast_column for statistical testing. Provide exactly 2 values.
 (str)
 (accepts 2 values)


### `--participant-label`, `--participant_label`

The label(s) of the participant(s) that should be analyzed. The label corresponds to sub-<participant_label> from the BIDS sec (so it does not include "sub-"). If this parameter is not provided all subjects should be analyzed. Multiple participants can be specified with a space separated list.
 (str)
 (accepts one or more values)


### `--exclude-participant-label`, `--exclude_participant_label`

The label(s) of the participant(s) that should be excluded. The label corresponds to sub-<participant_label> from the BIDS spec (so it does not include "sub-"). If this parameter is not provided all subjects should be analyzed. Multiple participants can be specified with a space separated list.
 (str)
 (accepts one or more values)


### `--derivatives`

Path(s) to a derivatives dataset, for folder(s) that contains multiple derivatives datasets
 (default: False)
 (accepts zero or more values)


### `--filter-spim`, `--filter_spim`

(default: suffix=SPIM extension=ome.zarr sample=brain)
 (accepts one or more values)


### `--filter-mri`, `--filter_mri`

(default: suffix=T2w extension=nii.gz datatype=anat)
 (accepts one or more values)


### `--wildcards-spim`, `--wildcards_spim`

(default: subject sample acquisition staining)
 (str)
 (accepts one or more values)


### `--wildcards-mri`, `--wildcards_mri`

(default: subject session acquisition run reconstruction suffix)
 (str)
 (accepts one or more values)


### `--path-spim`, `--path_spim`

 (str)


### `--path-mri`, `--path_mri`

 (str)


### `--pybidsdb-dir`, `--pybidsdb_dir`

Optional path to directory of SQLite databasefile for PyBIDS. If directory is passed and folder exists, indexing is skipped. If pybidsdb_reset is called, indexing will persist
 (Path)


### `--pybidsdb-reset`, `--pybidsdb_reset`

Reindex existing PyBIDS SQLite database
 (default: False)
 (flag)


### `--reset-db`, `--reset_db`

==SUPPRESS==
 (default: False)
 (flag)


### `--help-snakemake`, `--help_snakemake`

Options to Snakemake can also be passed directly at the command-line, use this to print Snakemake usage


### `--force-output`, `--force_output`

Force output in a new directory that already has contents
 (default: False)
 (flag)

