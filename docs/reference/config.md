# Configuration Options Reference

<!-- TODO: Add comprehensive configuration reference -->

Complete reference for SPIMquant configuration options.

## Configuration File Structure

```yaml
# TODO: Add complete configuration schema
```

## Input Configuration

<!-- TODO: Document input options -->

## Processing Parameters

<!-- TODO: Document processing parameters -->

## Output Configuration

<!-- TODO: Document output options -->

## Template Options

<!-- TODO: Document template configuration -->

## Segmentation Options

<!-- TODO: Document general segmentation configuration -->

### LANTERN plaque segmentation

Used when `--seg_method` includes `lantern`.  See [Segmentation Methods](../howto/segmentation.md#lantern-plaques-seg_method-lantern) for how these fit together.

| CLI option | Config key | Default | Meaning |
|---|---|---|---|
| `--lantern_config` | `lantern_config` | `ki3-abeta@v2` | Ensemble to run, as `<name>[@<revision>]` → Hugging Face repo `apooladi/lantern-<name>` |
| `--merge_mode` | `plaque_merge_mode` | `gaussian` | How overlapping tiles' votes are combined: `average`, `gaussian` or `max` |
| `--min_size` | `min_size` | `4` | Drop plaque components smaller than this many inference-grid voxels (0 disables) |
| `--label_filter` | `label_filter` | off | 3-class ensembles: drop plaques touching a large "false positive" component (UNVALIDATED) |
| | `stains_for_plaques` | `abeta`, `Abeta`, `BetaAmyloid` | LANTERN runs on the first of these present in the data |
| | `plaque_level` | `1` | Starting (and coarsest allowed) pyramid level for inference |
| | `plaque_iso_res` | `4.0` | Target isotropic inference resolution, µm |
| | `plaque_n_folds` | `5` | Folds in the ensemble |
| | `plaque_vote_threshold` | `3` | Folds that must vote plaque; with `average`/`gaussian`, a mean of at least this − 0.5 |
| | `plaque_fp_min_votes` | `1` | Folds that must call a voxel "false positive" (3-class ensembles) |
| | `label_filter_size` | `10` | Minimum false-positive component size, in voxels, for `--label_filter` |
| | `plaque_tile` / `plaque_stride` | `128` / `64` | Sliding-window tile size and step |
| | `plaque_batch_size`, `plaque_chunk`, `plaque_n_gpus` | `8`, `512`, `2` | Throughput settings: tiles per GPU batch, block size, GPUs per job |

## Resource Management

<!-- TODO: Document resource options -->

## Next Steps

- [Workflow Rules](rules.md)
- [Output Files](outputs.md)
