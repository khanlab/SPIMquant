# Group QC Report

The participant level writes QC figures and reports for each subject (in
`sub-*/qc/` and `sub-*/xfm/*_regqc.html`). The group level gathers them into one
page, `group/qc/groupqc.html`, where you can step through the subjects and
compare them with each other.

```bash
pixi run spimquant /path/to/bids /path/to/output group --cores all
```

The report is written with every group run (skip it with `--no_groupqc`). To
make only the report:

```bash
pixi run spimquant /path/to/bids /path/to/output group --until groupqc_report --cores all
```

Outputs:

| File | Contents |
|------|----------|
| `group/qc/groupqc.html` | The report: one HTML file to open in a browser |
| `group/qc/groupqc.tsv` | The metrics, one row per subject, with the flags |
| `sub-*/qc/*_groupqc.json` | The slices and metrics of each subject, which the report is built from |

## Views

Every image is a slice of a subject's image after registration to the template.
All subjects are shown at the same template slices, so differences between
subjects show up as changes in the image, with the atlas drawn as outlines.

- **Subjects**: one subject at a time. Its registration stain with atlas
  outlines, its other stains, its field fraction maps, and its metrics. Each
  metric has a strip plot of all subjects with this subject highlighted.
  Below these are the subject's own QC figures (registration report,
  segmentation slices, ROI montages, histograms, ...). Step through subjects
  with <kbd>←</kbd> <kbd>→</kbd> or the sidebar. Hold <kbd>T</kbd>, or the
  *Hold for template* button, to flick to the template at the same slices.
- **Montage**: every subject at one slice of one image (a stain or a field
  fraction map), ordered by name, number of flags, a `participants.tsv` column
  or a metric. Click a tile to open that subject.
- **Flip**: the same slice of each subject in turn, as an animation, with
  speed, pause and stepping (<kbd>Space</kbd>, <kbd>←</kbd> <kbd>→</kbd>), and
  optionally the template as the first frame. A poorly registered subject, or
  one with missing tissue or odd intensity, shows up as a jump.
- **Metrics**: the distribution of each metric and a sortable table of all
  subjects.

Field fraction maps can be colored on a scale shared by all subjects (the median
of the subjects' own maxima) or each on its own.

## Metrics and flags

| Metric | What it shows |
|--------|---------------|
| NCC to template | Correlation of the registered registration stain with the template image in the brain mask |
| NCC to group | Correlation with the mean of the other subjects at the report slices: low when a subject is placed or shaped differently from the rest |
| Brain coverage | Share of the template brain mask with image data: low when tissue is cut, missing or outside the field of view |
| Jacobian p1 / p99, folded voxels | Local compression and expansion of the warp, and where it folds |
| Median intensity, intensity p99 | Per stain, in the template brain mask (raw units) |
| Mean field fraction | Per segmented stain, from its template-space field fraction map |
| Objects, density, field fraction | Per segmented stain, over all atlas regions of the merged segstats table |

A metric is flagged when it is far from the rest of the group: a robust z-score
(median and MAD, with the spread at least 1% of the median) beyond
`groupqc: outlier_z` (3 by default) in the direction that signals a problem.
This needs at least 5 subjects. Flags point at what to look at; they do not
decide whether a subject is usable.

## Options

| Option | Effect |
|--------|--------|
| `--no_groupqc` | Do not write the report in group runs |
| `--groupqc_group_by COL [COL ...]` | `participants.tsv` columns to color, filter and order subjects by (default: the `--group-stats-pairwise` factors) |
| `--groupqc_embed_figures` | Embed the subjects' QC figures (downsized) so that the report can be shared on its own. By default they are linked, which keeps the report small but works only when it is opened in the output folder |

In the config, `groupqc:` sets the slices per plane, the image quality, the
flag threshold, whether to compute Jacobian metrics and the atlas drawn as
outlines (default `coarse`, else `roi22`). With the default 9 slices, the report
is about 300 KB per subject for two stains and one field fraction map, so about
15 MB for 50 subjects.
