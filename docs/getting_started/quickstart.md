# Quick start

This runs the small test dataset that comes with the repository: one mouse
brain, downsampled to about 100 µm, with Abeta, CD31 and YoPro channels. It
is too coarse for meaningful numbers, but it runs every step in about
<!-- QS-TIME -->20 minutes on 4 cores, and produces every kind of output.
[Example outputs](../examples/index.md) shows the same outputs from a
full-resolution brain.

## 1. Dry run

From the repository (see [Installation](installation.md)):

```bash
pixi run spimquant $PWD/tests/bids_ds /tmp/spimquant_out participant -n \
  --registration-level 0 --segmentation-level 0 --atlas-segs roi198
```

`-n` builds the list of jobs without running them. The test image has only a
few pyramid levels, so registration and segmentation both use level 0 here;
for real data leave the defaults. The start of the output says which stain
is used for registration and which method segments each stain:

```
Segmentation methods (first of each is primary):
  Abeta: seg_method gmm+n3k1 (default, no-GPU fallback), object_method components (default)
```

## 2. Run

```bash
pixi run spimquant $PWD/tests/bids_ds /tmp/spimquant_out participant --cores 4 \
  --registration-level 0 --segmentation-level 0 --atlas-segs roi198
```

## 3. Look at the results

```
/tmp/spimquant_out/
├── sub-AS36F2/
│   ├── micr/        # images at the registration level, brain mask, warped to the template
│   ├── xfm/         # warps, and the registration QC report (*_regqc.html)
│   ├── parc/        # atlas labels in the subject
│   ├── seg/         # masks (.ozx), field fraction and count maps
│   ├── vessels/     # CD31 vessel mask and distance map
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
pixi run spimquant $PWD/tests/bids_ds /tmp/spimquant_out group --cores 4 --atlas-segs roi198
```

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
