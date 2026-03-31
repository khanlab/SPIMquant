# Output Files Reference

<!-- TODO: Add comprehensive output files reference -->

Reference for SPIMquant output files and formats.

## Directory Structure

```
output/spimquant/
├── sub-01/
│   ├── micr/
│   │   ├── *_space-{template}_SPIM.nii.gz
│   │   ├── *_space-{template}_regqc.html
│   │   └── *_mask.nii.gz
│   ├── seg/
│   │   └── *_SPIM.ome.zarr
│   ├── parc/
│   │   └── *_from-{template}_dseg.nii.gz
│   ├── tabular/
│   │   ├── *_segstats.tsv
│   │   ├── *_regionpropstats.tsv
│   │   ├── *_fieldfracstats.tsv
│   │   └── *_countstats.tsv
│   └── xfm/
│       └── *_xfm.{txt,nii.gz}
└── group/
    ├── *_groupstats.tsv
    ├── *_groupstats.png
    └── *_groupstats.nii.gz
```

## Participant-Level Outputs

### Registered Images (`micr/`)

`*_space-{template}_SPIM.nii.gz`

Registered SPIM data in template space, stored as NIfTI files.

### Registration QC Reports (`micr/`)

`*_space-{template}_regqc.html`

HTML reports for visually checking registration quality.

### Brain Masks (`micr/`)

`*_mask.nii.gz`

Binary brain mask in subject space.

### Segmentation Data (`seg/`)

`*_SPIM.ome.zarr`

Full-resolution segmentation results stored as OME-Zarr arrays.

### Parcellation Maps (`parc/`)

`*_from-{template}_dseg.nii.gz`

Atlas-based discrete segmentation (parcellation) maps in subject space.

### Statistics Tables (`tabular/`)

`*_segstats.tsv`

Per-region segmentation statistics combining region properties, field fraction, and count metrics.

`*_regionpropstats.tsv`

Per-region object-level properties (e.g. volume, centroid).

`*_fieldfracstats.tsv`

Per-region field fraction (percentage of voxels positive) statistics.

`*_countstats.tsv`

Per-region object count statistics.

### Transforms (`xfm/`)

`*_xfm.txt` / `*_xfm.nii.gz`

Affine and deformable registration transforms between subject and template space.

## Group-Level Outputs

Group-level outputs are stored directly under `output/spimquant/group/` (no subject subdirectory).

### Statistical Results

`*_groupstats.tsv`

Statistical test results (t-statistics, p-values, effect sizes) for each brain region.

### Visualizations

`*_groupstats.png`

Heatmap visualizations of statistical results across brain regions.

### Volume Maps

`*_groupstats.nii.gz`

3D volumetric maps of statistical values for visualization in neuroimaging software.

## Quality Control Outputs

<!-- TODO: Document QC outputs -->

## Intermediate Files

<!-- TODO: Document intermediate files -->

## Next Steps

- [API Documentation](api.md)
- [Configuration Options](config.md)