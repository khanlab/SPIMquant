# Dataset QC Report

The `qc` analysis level checks a batch of scans before committing to
segmentation and quantification. It does no segmentation; it writes one
self-contained HTML report on the raw data, with a table of its metrics.

```bash
pixi run spimquant /path/to/bids /path/to/output qc --cores all
```

Outputs:

| File | Contents |
|------|----------|
| `qc/qcreport.html` | The report: one file with all images embedded, to open in a browser or send to someone |
| `qc/qcmetrics.tsv` | The metrics, one row per scan and stain |
| `sub-*/qc/*_qcscan.json` | The metrics and images of each scan, which the report is built from |

## Scans with different stains

Every other analysis level needs all scans to have the same stains. The `qc`
level does not: each scan is reported on with the stains it has, so that a
dataset piloting different stainings can be looked at as a whole. The report
has a column per stain, and a scan without that stain has an empty cell.

Stain labels that differ only in case share a column, as do labels listed
together under `qc: stain_aliases` in the config (e.g. `YoPro` and `YOPRO`).
`--ignore_stains` and `--sample_staining` work as elsewhere.

## What is in the report

- **Scans by stain**: a thumbnail and a metric of your choice for every scan
  and stain, with warnings on values that stand out. Select a cell to open it.
- **One scan**: overview slices in three planes with atlas outlines,
  full-resolution crops of every stain at the same places, the metrics, the
  intensity histogram and the intensity through the sample, and a pass / check /
  fail rating with a note.
- **One stain across scans**: the overview, or the crop of one region, of
  every scan that has the stain, side by side, optionally ordered by a column
  of `participants.tsv` such as a staining protocol (`--qc_group_by`).

All scans of a stain are shown with the same display window, which can be
adjusted, so brightness can be compared between scans.

Ratings are kept in the browser they are made in; download them as a TSV from
the report to keep or share them.

## Metrics

Per stain, from the image at `--registration_level`, within a tissue mask:

| Metric | Meaning |
|--------|---------|
| SNR | Tissue median above background, over the background noise |
| Tissue / background | Tissue median over background median |
| Signal range | 99.9th percentile over the median of tissue: high for sparse bright labelling, near 1 for flat signal |
| Depth non-uniformity | Variation of the median tissue intensity from slice to slice, along the worst axis (attenuation, uneven clearing) |
| Sharpness | Fine detail relative to coarse structure in the full-resolution crops |
| Saturated | Fraction of pixels at the top of the data type, in the worst crop |

Per scan: the tissue volume, the fraction of the tissue surface at the edge of
the field of view, and how well the registered atlas overlaps the tissue.

A value is flagged when it is far from the other scans of the same stain
(which needs at least four of them), or beyond a fixed limit. Flags say where
to look; they do not decide whether a scan is usable.

## Registration and crops

Each scan is registered to the template with a quick greedy registration of
its own first stain from `--stains_for_reg`, on the raw downsampled image.
This is separate from, and rougher than, the registration of the
`participant` level; it is only used to draw atlas outlines and to take the
crops in the same atlas regions in every scan. Add `--sloppy` to make it
faster still.

Scans without any of `--stains_for_reg` are not registered, and get crops
spread through the tissue instead. `--qc_no_reg` does this for all scans.

| Option | Default | |
|--------|---------|---|
| `--qc_atlas_seg` | `coarse`, if the template has it | Atlas the crop regions and outlines come from |
| `--qc_crop_labels` | the largest regions | Regions to crop, by name, abbreviation or index |
| `--qc_n_crops` | 6 | Number of crops per scan |
| `--qc_crop_size` | 512 | Size of the crops in pixels |
| `--qc_crop_level` | 0 | Pyramid level of the crops |

The size of the report grows with the number of scans, stains and crops;
around 12 scans of 4 stains with the defaults gives a file of a few tens of MB.
