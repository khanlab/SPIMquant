# Output Files Reference

<!-- TODO: Add comprehensive output files reference -->

Reference for SPIMquant output files and formats.

## Directory Structure

```
output/
└── spimquant/
    ├── sub-01/
    │   └── micr/
    │       ├── *_space-template_SPIM.nii.gz
    │       ├── *_dseg.nii.gz
    │       └── *_segstats.tsv
    └── qc/
```

## Participant-Level Outputs

<!-- TODO: Document participant outputs -->

### Registered Images

`*_space-template_SPIM.nii.gz`

<!-- TODO: Add file format details -->

### Segmentation Maps

`*_dseg.nii.gz`

<!-- TODO: Add file format details -->

### Statistics Tables

`*_segstats.tsv`

<!-- TODO: Add table format details -->

## Group-Level Outputs

<!-- TODO: Document group outputs -->

### Statistical Results

`*_groupstats.tsv`

<!-- TODO: Add file format details -->

### Visualizations

`*_groupstats.png`

<!-- TODO: Add visualization details -->

### Volume Maps

`*_groupstats.nii`

<!-- TODO: Add file format details -->

## Quality Control Outputs

<!-- TODO: Document QC outputs -->

## Intermediate Files

<!-- TODO: Document intermediate files -->

## Next Steps

- [API Documentation](api.md)
- [Configuration Options](config.md)