# Segmentation Methods

SPIMquant detects pathology signal (e.g., amyloid plaques, microglia) in high-resolution SPIM data
using intensity-based thresholding combined with optional k-means refinement.  This guide explains
each step, the available options, and how to tune them for your data.

For a general overview of where segmentation fits in the pipeline, see the
[Workflow Overview](../reference/workflow_overview.md#stage-4-pathology-segmentation-segmentationsmk).

---

## Step 1 — Bias-Field Correction

Before segmentation, low-frequency intensity gradients caused by the lightsheet illumination
profile are removed.  Two methods are available via `--correction_method`:

| Method | Flag | Speed | Quality |
|--------|------|-------|---------|
| Gaussian smoothing | `gaussian` | Fast | Good for most data |
| ANTsPy N4 | `n4` | Slower | Better for strong gradients |

```bash
# Use gaussian correction (fast)
pixi run spimquant /bids /out participant --correction_method gaussian

# Use N4 correction (more accurate, default)
pixi run spimquant /bids /out participant --correction_method n4
```

---

## Step 2 — Segmentation Algorithm

The `--seg_method` flag controls the segmentation approach:

### `threshold`

A global intensity threshold is applied to the bias-corrected volume.  This is the simplest
approach and works well when the pathology signal is well-separated from background.

### `otsu+k3i2` (default)

A multi-Otsu threshold is first computed from the image histogram.  The Otsu-derived threshold
separates candidate foreground voxels from background.  A k-means clustering step (3 classes,
2 iterations) is then applied to the candidate foreground region to further refine the
segmentation boundary.  This method is robust to overlapping intensity distributions.

The algorithm (in `multiotsu.py`):
1. Compute image histogram with preset bin range to handle intensity outliers
2. Use `skimage.filters.threshold_multiotsu` to find thresholds between *k* classes
3. Apply the selected threshold index to produce a binary mask
4. Multiply by 100 so that downsampling yields percent occupancy (field fraction)

---

## Step 3 — Object Filtering

Small objects and edge artefacts are removed by `compute_filtered_regionprops`:

- **Minimum area filter**: Objects smaller than `regionprop_filters.area_min` voxels
  (default: 200) are discarded.
- **Edge proximity filter**: Objects whose centroid is within a configurable distance of the
  brain mask edge are removed to eliminate illumination artefacts.

Configuration in `snakebids.yml`:

```yaml
regionprop_filters:
  area_min: 200   # minimum object size in voxels
```

---

## Output Mask Scaling

!!! warning "Mask values are 0 and 100, not 0 and 1"
    SPIMquant stores binary segmentation masks with values **0** (background) and **100**
    (foreground).  This is intentional: when the mask is downsampled using local-mean
    pooling to compute field fraction, the resulting values directly represent the percentage
    of voxels occupied by pathology signal (0–100%).

---

## Choosing Stains for Segmentation

SPIMquant automatically detects which stains to segment based on the `--stains_for_seg`
configuration (default: `abeta`, `Abeta`, `BetaAmyloid`, `AlphaSynuclein`, `Iba1`, `ChAT`).
Only stains present in your BIDS dataset and matching this list will be segmented.

```bash
# Override the default stain list
pixi run spimquant /bids /out participant \
  --stains_for_seg abeta Iba1
```

---

## Next Steps

- [Output Files Reference](../reference/outputs.md#stage-4-pathology-segmentation): All segmentation outputs
- [Imaris Crops](imaris_crops.md): Extract high-resolution region crops
- [Workflow Overview](../reference/workflow_overview.md): Full pipeline description