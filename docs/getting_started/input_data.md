# Input data

SPIMquant reads a BIDS dataset with one multi-channel lightsheet scan per
brain in the `micr` folder. [SPIMpack](https://github.com/khanlab/SPIMpack)
(formerly SPIMprep) writes this layout from raw Blaze, SmartSPIM or Imaris
data.

```
bids/
├── dataset_description.json
├── participants.tsv
├── samples.tsv
└── sub-01/
    └── ses-exvivo/                   # optional
        └── micr/
            ├── sub-01_ses-exvivo_sample-brain_acq-4x_SPIM.ims
            └── sub-01_ses-exvivo_sample-brain_acq-4x_SPIM.json
```

## Image formats

| Format | Extension | Selected by default |
|---|---|---|
| OME-Zarr (directory) | `.ome.zarr` | yes |
| Imaris (HDF5) | `.ims` | yes |
| OME-Zarr in a zip store | `.ozx`, `.ome.zarr.zip` | no: add `--filter-spim extension=ozx` (or `extension=ome.zarr.zip`) |

The image must hold a resolution pyramid: registration runs on a coarse
level (`--registration-level`, default 5) and segmentation on a fine one
(`--segmentation-level`, default 0, full resolution). Voxel sizes are read
from the image metadata. Imaris files are read through a process-safe reader,
so several jobs can read one file at once.

## Stains and orientation

Each channel needs a stain name. SPIMquant takes it, in order of precedence,
from `--sample_staining`, the `SampleStaining` list in the `_SPIM.json`
sidecar, or the OME-Zarr channel labels. The stain names decide what happens
to each channel:

| Option | Default stains | Used for |
|---|---|---|
| `--stains-for-reg` | PI, YOPRO, YoPro, YoPro1, AutoF, autof | registration (the first one present) |
| `--stains_for_seg` | abeta, Abeta, BetaAmyloid, AlphaSynuclein, Iba1, NeuN, Olig2, ChAT, VAChT, Thy1GFP | segmentation and per-region tables |
| `--stains_for_vessels` | CD31, Lectin | vessel segmentation (the first one present) |
| `--stains-for-fibres` | VAChT, vAChT, MOG, Mog | fibre orientation, with `--fibres` |

A stain named in `--seg_method STAIN=...` is segmented even if it is not in
`--stains_for_seg`. Use `--ignore-stains` to drop channels that only some
brains have, so that the cohort is uniform.

If the same dye is labelled differently in different datasets, map the
names with `--stain-synonyms ALIAS=CANONICAL` (or the `stain_synonyms` config
mapping). Channels labelled ALIAS are then processed, and named in the outputs,
as CANONICAL. For example, to register PI-stained brains to the YoPro cohort
template: `--stain-synonyms PI=YoPro --via-template YoPro99`.

The orientation of the image axes comes from `--orientation`, the
`OrientationStringXYZ` key of the sidecar, or the OME-Zarr metadata. It is
written in XYZ order with the letter each axis increases towards: Blaze data
stitched with SPIMpack is usually `RPI`, Imaris-stitched data `RAI`, and
LifeCanvas data `RAS`, but a brain mounted upside down changes this. A wrong
orientation shows up immediately in the registration QC report as a flipped
or rotated brain; the [dataset QC report](../howto/dataset_qc.md) lets you
check it for every scan before processing.

Example sidecar:

```json
{
  "OrientationStringXYZ": "RAI",
  "PixelSize": [0.001625, 0.001625, 0.00275],
  "PixelSizeUnits": "mm",
  "SampleStaining": ["Abeta", "CD31", "YoPro"]
}
```

## Participants and samples

`participants.tsv` holds one row per subject. Its columns (genotype,
treatment, age, sex...) are copied into the group tables and can be used in
the group model. `samples.tsv` lists the samples (`sample-brain`).

## Selecting scans

The usual snakebids options select what is processed:

```bash
--participant-label 01 02          # only these subjects
--exclude-participant-label 03
--filter-spim acquisition=4x       # only scans with acq-4x
```

A second scan of the same brain (e.g. a high-resolution scan of one region)
can be processed alongside the reference scan with `--aux_spim`; see
[Multiple scans of a sample](../howto/multiple_scans.md).
