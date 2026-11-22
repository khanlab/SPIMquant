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

### Segmentation methods per stain

See [Segmentation Methods](../howto/segmentation.md#choosing-methods-per-stain).

| CLI option | Config key | Default | Meaning |
|---|---|---|---|
| `--seg_method` | `seg_method` | none (built-in defaults) | Mask method(s) per stain: `STAIN=METHOD[,METHOD...]` or a bare `METHOD` for every stain not named; first is primary |
| `--object_method` | `object_method` | none (built-in defaults) | Object method(s) per stain, same syntax: `watershed`, `components` or `none` |
| | `seg_method_defaults` | `default: gmm+n3k1`; `abeta`/`Abeta`/`BetaAmyloid`: `lantern` | Built-in per-stain mask methods |
| | `object_method_defaults` | `default: watershed`; `abeta`/`Abeta`/`BetaAmyloid`: `components` | Built-in per-stain object methods |
| | `sdt_stains` | `[]` | Stains (besides the vessel stain) whose primary mask's signed distance is sampled at every object (`<stain>_sdt`) |

### LANTERN plaque segmentation

Used for every stain whose methods include `lantern` (by default amyloid-beta).  See [Segmentation Methods](../howto/segmentation.md#lantern-plaques-seg_method-lantern) for how these fit together.

| CLI option | Config key | Default | Meaning |
|---|---|---|---|
| `--lantern_config` | `lantern_config` | `abeta` | Ensemble to run, as `<name>[@<revision>]` → Hugging Face repo `apooladi/lantern-<name>`; `<revision>` defaults to `main` |
| `--merge_mode` | `plaque_merge_mode` | `gaussian` | How overlapping tiles' votes are combined: `average`, `gaussian` or `max` |
| `--min_size` | `min_size` | `4` | Drop plaque components smaller than this many inference-grid voxels (0 disables) |
| `--label_filter` | `label_filter` | off | 3-class ensembles: drop plaques touching a large "false positive" component (UNVALIDATED) |
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
