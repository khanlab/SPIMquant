# Quick start

This runs the small test dataset that comes with the repository: one mouse
brain, downsampled to about 100 µm, with Abeta, CD31 and YoPro channels. It
is too coarse for meaningful numbers, but it runs every step in about
35 minutes on 4 cores, and writes the same kinds of outputs as a real run.
[Example outputs](../examples/index.md) shows the same outputs from a
full-resolution brain.

## 1. Dry run

From the repository (see [Installation](installation.md)):

```bash
OPTS="--registration-level 0 --segmentation-level 0 --atlas-segs roi198"
pixi run spimquant $PWD/tests/bids_ds /tmp/spimquant_out participant -n $OPTS
```

`-n` builds the list of jobs without running them. The test image is tiny
(about 100 µm voxels, few pyramid levels), so these options run registration
and segmentation at level 0 and use a single atlas (`roi198`); for real data
leave them out. The start of the output says how the run is set up:

```
Note: using --correction-level 2 instead of 6: at that level an image has fewer than 16 voxels along an axis, too few for N4
Segmentation methods (first of each is primary):
  Abeta: seg_method gmm+n3k1 (default, no-GPU fallback), object_method components (default)
```

Without a GPU, Abeta is segmented with `gmm+n3k1` instead of LANTERN, and
vesselFM segments CD31 on the CPU (seconds on this image; a whole brain at
full resolution takes 2.9 h on one A100 GPU, and far longer on a CPU).

## 2. Run

```bash
pixi run spimquant $PWD/tests/bids_ds /tmp/spimquant_out participant --cores 4 $OPTS
```

## 3. Look at the results

```
/tmp/spimquant_out/
├── sub-AS36F2/
│   ├── micr/        # images at the registration level, brain mask, warped to the template
│   ├── xfm/         # warps, and the registration QC report (*_regqc.html)
│   ├── parc/        # atlas labels in the subject
│   ├── seg/         # masks (.ozx), field fraction and count maps
│   ├── vessels/     # CD31 vessel mask and field fraction
│   ├── tabular/     # per-region and per-object tables
│   ├── featuremap/  # per-region values painted into the template
│   └── qc/          # QC figures
├── tpl-ABAv3/       # the template and atlas used
└── logs/
```

Start with:

1. `sub-AS36F2/xfm/*_regqc.html`: open it in a browser; the atlas contours
   should follow the brain's anatomy.
2. `sub-AS36F2/qc/*_segslices.png`: the Abeta mask over the image.
3. `sub-AS36F2/tabular/*_seg-roi198_from-ABAv3_mergedsegstats.tsv`: one row
   per region, one column per stain and metric:

```python
import pandas as pd
df = pd.read_csv("sub-AS36F2/tabular/sub-AS36F2_sample-brain_acq-downsampled_seg-roi198_from-ABAv3_mergedsegstats.tsv", sep="\t")
df.sort_values("Abeta+fieldfrac", ascending=False)[["name", "Abeta+fieldfrac", "Abeta+density"]].head()
```

Every file is described in the [output reference](../reference/outputs.md).

## 4. Group level

```bash
pixi run spimquant $PWD/tests/bids_ds /tmp/spimquant_out group --cores 4 $OPTS
```

Give the group level the same options as the participant level: it checks
that the participant outputs are up to date, and with different options it
would plan (and run) a different participant level.

With one subject this only stacks the tables into
`group/1/seg-roi198_from-ABAv3_allsubjects.tsv`, with the columns of
`participants.tsv` added. With more subjects, add a model to compare groups;
see [Group analysis](../usage/group_analysis.md).

## Your own data

- Arrange it as described in [Input data](input_data.md), and check it with
  the [dataset QC report](../howto/dataset_qc.md).
- Drop the level options; run with `--cores all` on a workstation, or on a
  cluster as in [Running SPIMquant](../usage/running.md).
- For amyloid-beta, run in the `gpu` environment to use LANTERN.
