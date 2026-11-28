# Template-space data, crops and patches

The participant level quantifies in subject space and only warps coarse maps
to the template. These options also write full-resolution image data, either
resampled into the template or cut from atlas regions.

## Whole channels in template space (`--template-spim-res`)

Resamples whole channels into template space as OME-Zarr, at an isotropic
resolution in microns, e.g. close to the input's:

```bash
spimquant bids/ out/ participant --template-spim-res 2.5
```

The output grid covers the template's field of view with voxels of exactly
that size: with the 25 µm ABAv3 template, `2.5` gives a grid 10 times finer
(where the resolution does not divide the template's extent, the grid extends
past it by less than a voxel). Data are sampled from pyramid level
`--template-spim-level` (default 0) through the subject-to-template warp, and
written with pyramid levels down to about 128 voxels per axis.

Each subject, stain and resolution gives
`micr/*_stain-{stain}_space-{template}_res-2p5um_level-0_desc-raw_SPIM.ozx`
(`p` for a decimal point). With `--template_spim_corrected` the
intensity-corrected image is resampled instead (`desc-corrected<method>`).
`--template_spim_chunks` sets the output chunk size (default 256);
`--template_spim_io direct` lets every output chunk read its own input, which
may suit the distributed scheduler with object storage.

## Boxes in template space (`--bbox_config`)

Writes one NIfTI per box, subject and stain, resampled into template space.
The boxes are listed in a YAML file:

```yaml
bounding_boxes:
  - name: entorhinal      # used as the desc- of the output; unique
    x: 5.2                # box centre in template RAS coordinates
    y: -2.3
    z: -1.8
    units: mm             # or um
    size: [256, 256, 256] # voxels along R, A, S (default 256 each)
    resolution_um: 20     # isotropic output voxel size
    input_level: 0        # input pyramid level (default 0)
  - name: wholebrain
    use_brain_mask: true  # box from the template brain mask instead
    resolution_um: 100
    input_level: 2
```

```bash
spimquant bids/ out/ participant --bbox_config boxes.yml
```

## Patches and crops of atlas regions (`--save_patches`)

Samples random patches inside atlas regions, in subject space, from the raw
image, the corrected image and the segmentation mask, for training or
reviewing segmentation models:

| Option | Default |
|---|---|
| `--patch-size X Y Z` | 256 256 256 voxels |
| `--patch-atlas-segs` | `coarse` |
| `--patch_labels` | all labels (names, abbreviations or indices) |
| `--n-patches-per-label` | 5 |
| `--patch-seed` | 42 |
| `--no_patch_uint8` | patches are scaled to uint8 by robust percentiles unless set |

## Imaris crops of atlas regions

With `--save_patches`, the bounding box of each atlas region of
`--crop-atlas-segs` (default `coarse`) is also cut out of the
full-resolution scan as an Imaris file, for viewing a region at full
resolution in Imaris or syGlass. `--crop_labels` limits this to some regions
(names, abbreviations or indices).

The file names of all of these are in the [output reference](../reference/outputs.md).
