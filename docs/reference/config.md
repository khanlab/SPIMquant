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

### LANTERN plaque inference

When `seg_method` includes `lantern`, SPIMquant runs the 5-fold LANTERN
amyloid-beta plaque ensemble on the first available stain from
`stains_for_plaques`.

Relevant configuration keys:

- `plaque_level` — pyramid level at which LANTERN inference runs
- `plaque_tile` / `plaque_stride` — 3D sliding-window tile size and stride
- `plaque_vote_threshold` — downstream vote threshold used when binarizing the
  probability map
- `plaque_merge_mode` (CLI: `--merge_mode`) — how overlapping plaque tiles are
  merged:
  - `max` — historical/default behavior; take the maximum per-voxel fold-vote
    count over overlapping tiles
  - `average` — uniformly average per-tile fold-vote fractions
  - `gaussian` — Gaussian-weighted average that emphasizes tile centers over
    tile edges

See [Segmentation Methods](../howto/segmentation.md#lantern-plaques-seg_method-lantern)
for usage guidance and trade-offs between the merge modes.

## Resource Management

<!-- TODO: Document resource options -->

## Next Steps

- [Workflow Rules](rules.md)
- [Output Files](outputs.md)