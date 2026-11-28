# FAQ

### Do I need a GPU?

No. Without one, SPIMquant uses CPU methods everywhere: Abeta is segmented
with `gmm+n3k1` instead of LANTERN, and a note is printed. A GPU (and the
`gpu` pixi environment) is needed for LANTERN, vesselFM vessel segmentation,
FireANTs registration and building cohort templates. See
[Running SPIMquant](usage/running.md#gpu-steps).

### Which level should I segment at?

Level 0 (full resolution, the default) for objects such as plaques and cells.
Higher levels are coarser and much faster: each level halves the resolution,
at least in-plane, so level 1 has at most a quarter of the voxels of level 0.
Object counts depend on the level, because small objects merge or vanish at
coarser levels, so compare brains segmented at the same level.

### Registration failed or looks wrong. What do I check first?

Open `xfm/*_regqc.html`. A brain that is flipped or rotated by 90° means the
orientation is wrong (`--orientation`, or `OrientationStringXYZ` in the
sidecar; see [Input data](getting_started/input_data.md#stains-and-orientation)).
If only part of the brain was imaged, crop the template to one hemisphere with
`--template-crop left|right`. If the contrast differs a lot from the
template's, try another `--registration_method` or a
[cohort template](howto/cohort_template.md).

### A stain is not segmented.

Check that its name is in `--stains_for_seg` (case matters), or name it in
`--seg_method STAIN=METHOD`. The stains and methods are printed at the start
of every run.

### Can different stains use different methods?

Yes: `--seg_method Abeta=lantern Iba1=otsu+k3i2 --object_method Iba1=watershed`.
See [Segmentation](howto/segmentation.md#choosing-methods-per-stain).

### Which file do I use for statistics?

Per subject, `tabular/*_mergedsegstats.tsv` (one row per region, a column per
stain and metric). For the cohort, run the `group` level and use
`group/<label>/*_allsubjects.tsv`, which adds the columns of
`participants.tsv`. Object-level analyses start from
`tabular/*_regionprops.parquet`.

### Are the group p-values corrected for multiple comparisons?

No. See [Group analysis](usage/group_analysis.md#region-wise-statistics).

### Can I reuse registrations or segmentations from an earlier run?

Yes, with `--reuse-derivatives`. Brain masks and transforms are reused by
default; segmentation and vessel masks with `--reuse_outputs ... segmentation
vessels`. See [Reusing outputs](howto/reuse_derivatives.md).

### "Lock directory exists"

A previous run was interrupted. If no other run is using the output
directory, add `--unlock` once, then run again.

### How do I report a problem?

Open an issue at <https://github.com/khanlab/SPIMquant/issues> with the
command (`spimquant --history OUT` shows it), the version
(`spimquant --version`) and the log of the failed job (`OUT/logs/`).
