# Adding a Deep-Learning Model

SPIMquant runs deep-learning segmentation models (LANTERN for amyloid-beta plaques, vesselFM for vessels) through one shared tiled-inference engine, `spimquant/workflow/scripts/tiled_inference.py`.  A new model doesn't need its own tiling, parallelisation or I/O code.  It needs a small **adapter** that says how to load the model and how to turn a batch of image tiles into predictions, plus a Snakemake rule.

This guide covers what the engine does, what an adapter must provide, the templates to start from, and how to test, tune and integrate the result.

## What the engine does

For a whole-brain image, `run_inference` does the following:

1. **Plans blocks.**  The image is cut into `chunk`-sized blocks (512³ by default).  Each block is read with a halo of `tile - stride` voxels, so tiles overlap correctly across block boundaries and there are no seams at block edges.  Blocks are independent: each writes only its own core.
2. **Skips air.**  A block that doesn't touch the low-resolution brain mask is never read, and a tile with no brain is never run.  Most of a light-sheet frame is background, so this is usually the biggest saving.
3. **Reads in parallel.**  Blocks are read by a pool of separate processes (`inference_read_procs`, 8 by default).  HDF5, which backs Imaris files, allows only one read at a time per process, so threads, including dask's threaded scheduler, can't read Imaris in parallel.  Separate processes can, and reads over NFS scale with their number.
4. **Tiles, batches and merges on the GPU.**  Each block is copied to the GPU once.  `tile`³ tiles at `stride` are run in batches of `batch_size`, and overlapping tiles are merged on the GPU (`merge_mode`: `gaussian`, `average` or `max`).
5. **Uses every GPU.**  The model is loaded once per GPU (`n_gpus`).  A block borrows a whole GPU while it runs, so GPUs share the work without a fixed split.
6. **Hands back the result.**  The merged output, passed through the adapter's `finalize`, is returned as a lazy dask array for the script to write as OME-Zarr.

## The adapter

An adapter subclasses `tiled_inference.ModelAdapter`:

| Member | Required | What it does |
|---|---|---|
| `out_labels` | yes | One label per output channel.  The output gets one channel each. |
| `load(device)` | yes | Load the model onto `device` and return it (or anything `predict_tiles` needs).  Called once per GPU, before any block. |
| `predict_tiles(model, tiles)` | yes | `tiles` is a `(B, t, t, t)` float32 batch of **raw** intensities, already on the model's device.  Return `(n_out, B, t, t, t)` per-tile outputs, ideally as a torch tensor on that device.  Preprocessing (normalisation) happens here. |
| `tile_dtype` | no (float32) | dtype of `predict_tiles`' output.  It's kept as-is under `merge_mode="max"`; the averaging modes always merge in float32. |
| `finalize(merged)` | no (float32 as-is) | Turn the merged `(n_out, z, y, x)` dask array into what gets stored, e.g. threshold to a mask or divide vote counts by the fold count. |
| `describe()` | no | One line for the run log. |

The two adapters already in the repository are worked examples:

- **`VesselFM`** in `vesselfm.py` is a single network.  It rescales each tile's 1st–99th percentile, outputs a sigmoid probability, and thresholds it after merging into a 0/100 mask.
- **`LanternEnsemble`** in `lantern_plaques.py` is a multi-fold ensemble with any number of folds, one shared encoder pass per batch, and per-tile conditioning statistics.  Its `predict_tiles` returns uint8 vote counts and `finalize` divides by the number of folds.

## Step by step

### 1. Write the inference script from the template

Copy [`docs/howto/templates/tiled_model.py`](templates/tiled_model.py) to `spimquant/workflow/scripts/<yourmodel>_seg.py` and fill in the `TODO`s:

- `build_network()`: construct the network.  Import it from the model's package, or vendor the definition as `lantern_model.py` does.
- `load()`: adapt to the checkpoint layout, e.g. `ckpt["state_dict"]`.
- `predict_tiles()`: apply the normalisation the model was **trained with**, run the network, and pick the output channel(s).  Return probabilities: they average sensibly across overlapping tiles, while logits and hard labels don't.
- `finalize()`: produce a 0/100 uint8 mask, or a float32 probability map if a later rule should do the thresholding (as LANTERN does, so the threshold can change without re-running inference).
- `main()`: choose how the input is read.  Either a pyramid level made near-isotropic (the template's default, as vesselFM), or a physical resolution the model was trained at (`target_scale` in µm, as LANTERN).

!!! warning "Don't name the script after a package it imports"
    A script called `vesselfm.py` that does `from vesselfm.zarrnii_plugin import ...` shadows the `vesselfm` package whenever the scripts folder is first on `sys.path`, e.g. in tests and notebooks.  Use a distinct name such as `<yourmodel>_seg.py`, and in tests append the scripts folder to `sys.path` rather than inserting it first.

### 2. Add the rules and config

Copy [`docs/howto/templates/tiled_model.smk`](templates/tiled_model.smk) into a rules file and replace `mymodel` throughout.  It has two rules:

- `import_mymodel_model` downloads the weights once into `resources/models/`, so inference jobs never touch the network.
- `run_mymodel` runs the script.  Keep `tasks_per_gpu=0`: without it the SLURM executor asks for one task per GPU, and `srun` launches the whole job once per GPU.

Add the weights URL and the tiling options to `spimquant/config/snakebids.yml`:

```yaml
models:
  mymodel: https://huggingface.co/<org>/<repo>/resolve/main/<file>.pt

mymodel_tile: 128        # the patch size the model was trained on
mymodel_stride: 64       # overlap is tile - stride; also the halo per block
mymodel_batch_size: 8
mymodel_chunk: 512       # block edge; memory per in-flight block grows with it
mymodel_merge_mode: gaussian
mymodel_threshold: 0.5
mymodel_n_gpus: 1
```

`inference_read_procs` is shared by all models.  Add `--<option>` entries to the `parse_args` section of the same file for the options users should be able to set on the command line.

### 3. Hook it into the pipeline

There are two ways to do this:

- **A separate output**, like vesselFM's `vessels/..._desc-vesselfm_mask.ozx`.  Add the output to the targets in the `Snakefile`, plus any downstream rules that consume it.
- **A segmentation method**, selectable with `--seg_method STAIN=<name>`, like LANTERN.  Downstream rules (field fraction, objects, QC, ...) request `seg/..._stain-<stain>_level-<segmentation_level>_desc-<name>_mask.ozx`, so a rule has to produce that file.  You also need to:
    - add `<name>` to `SEG_METHOD_PATTERNS` in `spimquant/workflow/scripts/seg_plan.py`;
    - add `<name>: <cpu fallback>` to `gpu_seg_method_fallbacks` in the config, so runs without a GPU fall back to a CPU method instead of failing;
    - optionally, make it a stain's default under `seg_method_defaults`;
    - restrict the rule's `stain` wildcard to the stains the model was trained on.  This matches how the LANTERN rules handle it.

    LANTERN runs on its own near-isotropic grid, writes a probability map, and has `binarize_lantern_plaques` threshold it and resample it to `segmentation_level`.  A model that runs at `segmentation_level` and thresholds in `finalize` can write the mask directly.

### 4. Test

- **Unit tests.**  `tests/test_tiled_model_template.py` shows the pattern.  It fills `build_network` with a tiny network, saves a checkpoint, and runs the adapter through `predict_tiles`, `finalize` and `run_inference` on a small OME-Zarr, on the CPU.  Also test your preprocessing against the reference implementation, as `test_vesselfm_scaling_matches_monai` does for vesselFM.
- **Against the reference implementation on a crop.**  Before a whole-brain run, run the model's own inference code and your adapter (`tiled_inference.predict_block` on the crop plus its halo) on the same crop of real data.  Compare the foreground fraction and Dice, and look at a slice.  This catches normalisation and orientation mistakes in minutes.
- **A whole subject.**  Run up to your rule on a GPU node:

    ```console
    pixi run -e gpu spimquant /bids /output participant --gpu yes --cores 32 \
        --resources gpu=2 --until run_mymodel
    ```

    The log prints the block plan (`N of M blocks touch brain`), progress with an ETA per block, and the tiles run versus skipped.

## Tuning

| Option | Effect |
|---|---|
| `tile` | Fixed by the model: use the patch size it was trained on. |
| `stride` | Overlap is `tile - stride`.  Each voxel is covered by `(tile / stride)³` tiles: `stride = tile/2` costs about 8× the GPU time of no overlap, and removes the tile-edge artifacts. |
| `merge_mode` | `gaussian` (weights tile centres most; the default), `average` or `max`. |
| `batch_size` | Larger batches use the GPU better, until GPU memory runs out. |
| `chunk` | Block edge.  Larger blocks re-read less halo, but use more memory per block in flight. |
| `n_gpus` | Blocks spread across GPUs; inference time divides roughly by the GPU count. |
| `inference_read_procs` | More readers help when reads are the bottleneck (low GPU use in `nvidia-smi`).  Each one costs memory (below). |

**Which is the bottleneck?**  Watch `nvidia-smi` during a run.  Near-constant 100% means the GPU is the limit: add GPUs or increase `stride`.  Frequent idle gaps mean the GPU is waiting on reads: add readers.

**Sizing resources.**  Some measurements from a 4x Imaris subject (A100 GPUs, 512³ chunks, 64-voxel halo), to start from and then refit for your model:

- each block in flight holds about 3 GB of host memory (its input plus float32 merge buffers), and about `(n_out + 2) × 4` bytes per voxel on the GPU (about 2 GB at 640³);
- each reader process holds about 1.6 GB (its own copy of the input's dask graph);
- vesselFM at level 0: about 17 s per block on one GPU and 8 s on two, plus about 27 minutes to write the pyramid, which doesn't speed up with more GPUs;
- Imaris reads over NFS: about 11 MB/s from 1 reader, 20 from 4 and 47 from 16.

## Pitfalls

- **Wrong normalisation or scale.**  Models are trained at a particular resolution and intensity transform.  Match both.  Read raw intensities unless the model was trained on corrected data (LANTERN and vesselFM both read raw).
- **Per-tile normalisation in air.**  Normalising each tile on its own makes pure-background tiles look like tissue (noise stretched to unit variance), which can produce false positives far from the brain.  The brain mask avoids running those tiles; keep it unless the model handles background.
- **Hard labels before merging.**  Merge probabilities, then threshold in `finalize`.  Thresholding each tile first loses the benefit of the overlap.
- **Reduced precision.**  `torch.autocast` with bfloat16 is faster on A100s, but check its output against fp32 on a crop before keeping it.
- **Using one model from several blocks at once.**  Don't share one model across GPUs or threads yourself.  `DevicePool` gives each block a whole GPU and its own copy of the model.
