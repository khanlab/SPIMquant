"""5-fold LANTERN ensemble Abeta plaque segmentation.

Runs a published LANTERN ensemble (``apooladi/lantern-<name>``, chosen with
``--lantern_config``; ``ki3-abeta`` by default) over the Abeta channel and
writes the per-fold agreement map -- the fraction of folds voting foreground,
so ``{0, .2, .4, .6, .8, 1}`` for five folds -- as an OME-Zarr.  Binarization
is deliberately left to a downstream rule so the vote threshold can be changed
without re-running inference.

Every fold is a decoder fine-tuned on the same frozen self-supervised encoder,
so the encoder is run once per batch and its skip features are fed to each
fold's decoder.  The encoder is ~2/3 of the FLOPs of a full forward pass, so
this roughly halves inference time against five independent passes, and it is
exact: the encoder is deterministic in eval mode and the decoders only read the
skips.  ``load_folds`` checks the encoders really are identical rather than
trusting the model card.

Deviations from ``vesselfm.py``, and why:

- ``ZarrNii.segment()`` cannot be used.  It dispatches through ``da.map_blocks``
  with disjoint blocks, no block position and a hardcoded ``uint8`` output, so it
  can express neither the stride-64 tile overlap the model requires nor
  brain-mask-based block skipping.  ``da.map_overlap`` is driven directly here.

- The model reads RAW intensities.  It was never trained on N4-corrected or
  otherwise rescaled data, so this path takes ``input.spim`` directly rather than
  the ``desc-corrected*`` store the GMM/Otsu rules consume.

- Blocks are driven by an explicit loop (``run_blocks``) rather than
  ``da.map_overlap``.  With map_overlap every input block is read before the
  brain test can run, so the air that makes up most of a frame is still pulled
  through the (slow, single-stream) Imaris/HDF5 reader.  The loop consults the
  low-res brain mask first and never reads a block that has no brain in it;
  the blocks it does read are fetched concurrently from a thread pool so that
  I/O overlaps GPU inference, and their votes land in a temporary zarr that
  the pyramid writer then reads back.

- The model is scale-sensitive: it was trained at ~4 um isotropic.  The input
  is loaded via ``load_near_isotropic`` rather than at a fixed pyramid level,
  because a fixed level lands on the right grid only for pyramids that
  downsample x/y alone (ki3); a pyramid that also halves z per level (mapt)
  puts z at 8 um by level 1, so the loader walks to a finer level and
  downsamples per axis to get near ``plaque_iso_res``.

Inference contract (from the model cards, must match training exactly):
  tiles of 128^3 at stride 64, per-tile 0.5/99.5 percentile clip then z-score,
  softmax over the classes, plaque where class-1 probability >= 0.5 per fold,
  overlapping tiles combined by MAX over votes.  3-class ensembles (classes
  background / plaque / "annotated false positive") additionally emit, as a
  second output channel, the number of folds whose argmax is class 2 -- the
  model card's rule is "class 2 where any fold's argmax is class 2", which is
  that channel at >= 1 vote, and binarize_lantern_plaques applies it.
"""

import contextlib
import itertools
import multiprocessing
import os
import queue
import tempfile
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import dask.array as da
import numpy as np
import torch
import zarr
from scipy.ndimage import binary_dilation
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii

from dask_setup import get_dask_client
from lantern_model import ResEncLUNet
from zarrnii_compat import (
    axis_scale_um,
    load_near_isotropic,
    open_inference_array,
    with_channels,
)


def tile_origins(extent, tile, stride):
    """Tile origins covering [0, extent), with the last pulled back to hit the edge."""
    if extent <= tile:
        return [0]
    origins = list(range(0, extent - tile + 1, stride))
    if origins[-1] + tile < extent:
        origins.append(extent - tile)
    return origins


def normalize_tiles(x):
    """Per-tile percentile clip then z-score, matching LANTERN's training transform."""
    out = torch.empty_like(x)
    quantiles = torch.tensor([0.005, 0.995], device=x.device, dtype=x.dtype)
    for i in range(x.shape[0]):
        lo, hi = torch.quantile(x[i].reshape(-1), quantiles)
        clipped = torch.clamp(x[i], lo, hi)
        std = clipped.std()
        out[i] = (clipped - clipped.mean()) / (std if std > 1e-8 else 1.0)
    return out


class TileCounter:
    """Thread-safe tally of tiles run vs. skipped, for the end-of-run summary."""

    def __init__(self):
        self.run = 0
        self.skipped = 0
        self._lock = threading.Lock()

    def add(self, run, skipped):
        with self._lock:
            self.run += run
            self.skipped += skipped


def predict_volume(
    vol, nets, device, tile, stride, batch_size, brain=None, counter=None
):
    """Tiled 5-fold vote maps for one 3D block, as uint8 in [0, len(nets)].

    Returns an array of shape (n_out, *vol.shape): channel 0 counts the folds
    that call plaque (class-1 probability >= 0.5); for a 3-class ensemble
    channel 1 counts the folds whose argmax is class 2.

    `brain`, if given, is a boolean array on the same grid as `vol`; tiles
    that contain no brain voxel are skipped and vote zero.  This is both the
    main throughput lever (most of a light-sheet frame is air) and a
    correctness fix: the per-tile clip + z-score normalisation amplifies the
    sensor noise of a pure-air tile to unit variance, which is out of
    distribution for a model trained on tissue patches and yields spurious
    votes far from the brain.

    The caller owns `device` exclusively for the duration of this call (see
    DevicePool), so no additional locking is needed here.
    """
    orig_shape = vol.shape
    pad = [max(0, tile - n) for n in orig_shape]
    if any(pad):
        # Blocks at the volume border can be thinner than one tile; pad, then crop back.
        vol = np.pad(vol, [(0, p) for p in pad], mode="edge")
        if brain is not None:
            brain = np.pad(brain, [(0, p) for p in pad], mode="edge")

    n_out = 2 if nets[0].n_classes >= 3 else 1
    votes = np.zeros((n_out,) + vol.shape, dtype=np.uint8)
    coords = [
        (z, y, x)
        for z in tile_origins(vol.shape[0], tile, stride)
        for y in tile_origins(vol.shape[1], tile, stride)
        for x in tile_origins(vol.shape[2], tile, stride)
    ]
    n_all = len(coords)
    if brain is not None:
        coords = [
            (z, y, x)
            for z, y, x in coords
            if brain[z : z + tile, y : y + tile, x : x + tile].any()
        ]
    if counter is not None:
        counter.add(len(coords), n_all - len(coords))
    if not coords:
        return votes[:, : orig_shape[0], : orig_shape[1], : orig_shape[2]]

    for start in range(0, len(coords), batch_size):
        batch = coords[start : start + batch_size]
        arr = np.stack(
            [vol[z : z + tile, y : y + tile, x : x + tile] for z, y, x in batch]
        )
        tensor = torch.from_numpy(arr).to(device)
        tensor = normalize_tiles(tensor)[:, None]
        batch_votes = torch.zeros(
            (n_out, len(batch), tile, tile, tile), dtype=torch.uint8, device=device
        )
        with (
            torch.no_grad(),
            torch.autocast(
                device.type, dtype=torch.bfloat16, enabled=device.type == "cuda"
            ),
        ):
            # One encoder pass serves every fold: load_folds has verified the
            # encoders are identical and pointed all folds at the same module.
            skips = nets[0].encoder(tensor)
            for net in nets:
                prob = torch.softmax(net.decoder(skips).float(), dim=1)
                batch_votes[0] += (prob[:, 1] >= 0.5).to(torch.uint8)
                if n_out == 2:
                    batch_votes[1] += (prob.argmax(dim=1) == 2).to(torch.uint8)
        batch_votes = batch_votes.cpu().numpy()

        for j, (z, y, x) in enumerate(batch):
            region = votes[:, z : z + tile, y : y + tile, x : x + tile]
            # MAX, not sum: a plaque clipped at one tile edge is recovered by the
            # tile that contains it whole.
            #
            # Folds are summed WITHIN a tile (above) and tiles are combined with max
            # (here). Those two steps do not commute -- max-then-sum would score a
            # voxel higher when different folds find it from different tiles -- but
            # this order is the one LANTERN's own inference uses
            # (scripts/infer_wholebrain.py), and plaque_vote_threshold was calibrated
            # under it. Swapping them would silently change what "3 of 5" means.
            np.maximum(region, batch_votes[:, j], out=region)

    return votes[:, : orig_shape[0], : orig_shape[1], : orig_shape[2]]


def share_encoder(nets, model_paths):
    """Point every fold at fold 0's encoder module, after checking they match.

    LANTERN fine-tunes only the decoder, on an encoder frozen from
    self-supervised pretraining, so every fold of an ensemble carries the same
    encoder weights.  That is what lets ``predict_volume`` run the encoder once
    per batch.  It is checked tensor-for-tensor rather than assumed: an ensemble
    trained with the encoder unfrozen would otherwise silently decode the wrong
    features.  Sharing the module also drops the redundant copies from device
    memory (the encoder is ~0.7 GB of fp32 weights per fold).
    """
    shared = nets[0].encoder
    reference = shared.state_dict()
    for net, path in zip(nets[1:], model_paths[1:]):
        state = net.encoder.state_dict()
        same = state.keys() == reference.keys() and all(
            torch.equal(state[k], reference[k]) for k in reference
        )
        if not same:
            raise ValueError(
                f"{path}: encoder weights differ from {model_paths[0]}; the folds "
                "of a LANTERN ensemble must share one frozen encoder"
            )
        net.encoder = shared
        # UNetDecoder keeps its own reference to the encoder (used only for
        # feature-map size bookkeeping); repoint it so the old copy is freed.
        net.decoder.encoder = shared


def load_folds(model_paths, device):
    """Instantiate every fold on `device`, once, in the main process.

    With the encoder shared (see ``share_encoder``) this is ~1.2 GB of fp32
    weights for a 5-fold ensemble, so replicating it per GPU is cheap.
    DevicePool then guarantees only one dask block uses a given device at a
    time, which bounds peak VRAM without needing a lock inside the forward pass.

    ``weights_only=True`` because these checkpoints are downloaded over the
    network; they carry only tensors and plain dicts, so nothing is lost.
    """
    nets = []
    for path in model_paths:
        ckpt = torch.load(path, map_location="cpu", weights_only=True)
        if ckpt.get("model_type") not in (None, "ResEncLUNet"):
            raise ValueError(
                f"{path} declares model_type={ckpt['model_type']!r}, "
                "but this script only builds ResEncLUNet"
            )
        n_classes = checkpoint_num_classes(ckpt, path)
        net = ResEncLUNet(num_classes=n_classes, num_input_channels=1)
        net.load_state_dict(ckpt["state_dict"])
        net.n_classes = n_classes
        nets.append(net.eval())
    if len({net.n_classes for net in nets}) != 1:
        raise ValueError(
            f"folds disagree on the number of classes: "
            f"{[net.n_classes for net in nets]}"
        )
    share_encoder(nets, list(model_paths))
    return [net.to(device) for net in nets]


def checkpoint_num_classes(ckpt, path):
    """Number of output classes of a checkpoint: 2 (binary) or 3.

    Taken from the checkpoint's own ``num_classes`` entry, cross-checked
    against the segmentation head's weight shape; older binary checkpoints
    carry no entry and are read from the head alone.
    """
    heads = [
        v.shape[0]
        for k, v in ckpt["state_dict"].items()
        if k.startswith("decoder.seg_layers.") and k.endswith(".weight")
    ]
    if not heads or len(set(heads)) != 1:
        raise ValueError(f"{path}: cannot read the class count from the seg head")
    from_head = int(heads[0])
    declared = ckpt.get("num_classes")
    if declared is not None and int(declared) != from_head:
        raise ValueError(
            f"{path} declares num_classes={declared} but its seg head has "
            f"{from_head} outputs"
        )
    if from_head not in (2, 3):
        raise ValueError(f"{path}: {from_head} classes; expected 2 or 3")
    return from_head


def resolve_devices(requested):
    """The devices to run on, clamped to what is actually visible."""
    if not torch.cuda.is_available():
        print("WARNING: no CUDA device visible; falling back to CPU", flush=True)
        return [torch.device("cpu")]
    available = torch.cuda.device_count()
    if requested > available:
        print(
            f"WARNING: {requested} GPUs requested but only {available} visible; "
            f"using {available}",
            flush=True,
        )
    return [torch.device(f"cuda:{i}") for i in range(max(1, min(requested, available)))]


class DevicePool:
    """Lends out one GPU at a time, with the full ensemble resident on each.

    The ensemble is replicated per device (~1.2 GB of fp32 weights for 5 folds
    sharing one encoder), so each GPU works independently on a whole dask block.  Handing out devices from
    a queue rather than striping blocks by index balances the load for free:
    blocks vary enormously in cost because most of a frame is background, so a
    static assignment would leave one card idle -- the same imbalance LANTERN's
    own two-GPU split had to correct for by splitting on cumulative tile count.
    """

    def __init__(self, devices, model_paths):
        self._free = queue.Queue()
        self.nets = {}
        for device in devices:
            self.nets[device] = load_folds(model_paths, device)
            self._free.put(device)

    @contextlib.contextmanager
    def acquire(self):
        device = self._free.get()
        try:
            yield device, self.nets[device]
        finally:
            self._free.put(device)


def make_brain_sampler(mask_path, image_shape):
    """A function mapping (z, y, x) slices on the inference grid to a brain mask.

    The mask lives at a much coarser pyramid level, so it is held in memory and
    nearest-neighbour sampled into whatever region is asked for, from the
    region's array location.  It is first dilated by one low-res voxel (tens of
    microns on the inference grid) so that a tile at the pial surface is never
    skipped on account of a slightly tight mask.
    """
    mask_znimg = ZarrNii.from_nifti(mask_path, axes_order="ZYX")
    mask_arr = np.asarray(mask_znimg.data).squeeze() > 0
    if mask_arr.ndim != 3:
        raise ValueError(f"expected a 3D brain mask, got shape {mask_arr.shape}")
    mask_arr = binary_dilation(mask_arr, iterations=1)

    scales = [m / i for m, i in zip(mask_arr.shape, image_shape[-3:])]

    def brain_of(slices):
        index = []
        for axis, sl in enumerate(slices):
            # centre of each hires voxel, mapped into the low-res grid
            src = np.floor((np.arange(sl.start, sl.stop) + 0.5) * scales[axis])
            index.append(np.clip(src.astype(int), 0, mask_arr.shape[axis] - 1))
        return mask_arr[np.ix_(*index)]

    return brain_of


def plan_chunks(shape, block, halo):
    """Chunk sizes for `shape` at `block`, with no chunk shorter than `halo`.

    A plain rechunk leaves a remainder chunk that can be far shorter than the
    halo (z=1030 at block 512 leaves a 6-voxel tail), and dask refuses to overlap
    a chunk smaller than the depth.  A short tail is absorbed into its
    predecessor instead.  Every chunk start stays a multiple of `block`, so the
    per-block tile grid stays aligned with the global one.
    """
    chunks = [(shape[0],)]
    for extent in shape[1:]:
        if extent <= block:
            chunks.append((extent,))
            continue
        full, remainder = divmod(extent, block)
        if remainder == 0:
            chunks.append((block,) * full)
        elif remainder < halo:
            chunks.append((block,) * (full - 1) + (block + remainder,))
        else:
            chunks.append((block,) * full + (remainder,))
    return tuple(chunks)


def block_plan(shape, chunks, halo):
    """(core, outer, inner) slice triples for every spatial block.

    `core` is the block on the full grid, `outer` is core grown by `halo` on
    each side (clipped to the volume), and `inner` is core expressed relative
    to outer, i.e. what to crop out of a prediction made on outer.  The halo is
    `tile - stride`, so every voxel of the core is covered by at least one tile
    lying wholly inside outer -- which is what makes the stride-64 overlap
    correct across block boundaries and not just within them.  An axis that
    is a single chunk gets no halo: there is no block boundary to bridge.
    """
    starts = [np.cumsum((0,) + tuple(axis[:-1])) for axis in chunks]
    plan = []
    for idx in itertools.product(*(range(len(axis)) for axis in chunks)):
        core, outer, inner = [], [], []
        for axis, i in enumerate(idx):
            lo = int(starts[axis][i])
            hi = lo + int(chunks[axis][i])
            h = 0 if len(chunks[axis]) == 1 else halo
            olo, ohi = max(0, lo - h), min(shape[axis], hi + h)
            core.append(slice(lo, hi))
            outer.append(slice(olo, ohi))
            inner.append(slice(lo - olo, hi - olo))
        plan.append((tuple(core), tuple(outer), tuple(inner)))
    return plan


_READER = {}


def _reader_init(spec, n_threads):
    """Open the inference array once in a reader process (see BlockReader)."""
    import dask

    dask.config.set(scheduler="threads", num_workers=n_threads)
    _READER["data"] = open_inference_array(**spec)


def _read_block(outer):
    data = _READER["data"]
    return np.asarray(data[(slice(None),) + outer])


class BlockReader:
    """Reads (c, z, y, x) blocks in worker processes.

    HDF5 (Imaris input) serialises readers behind one process-wide lock, so
    threads in a single process cannot exceed one stream (~20-40 MB/s on the
    cluster's NFS); separate processes each get their own stream and scale
    to what the storage delivers.  Each worker opens the array once, at
    start, from `spec` (the arguments to open_inference_array); a request
    then ships only the slice tuple out and the block back.  Shipping the
    dask graph itself per block, as a distributed client would, costs
    ~100 MB of pickling per request on an Imaris pyramid and starves the
    workers.
    """

    def __init__(self, spec, n_procs, threads_per_proc=2):
        self._pool = ProcessPoolExecutor(
            max_workers=n_procs,
            mp_context=multiprocessing.get_context("spawn"),
            initializer=_reader_init,
            initargs=(spec, threads_per_proc),
        )

    def __call__(self, outer):
        return self._pool.submit(_read_block, outer).result()

    def close(self):
        self._pool.shutdown(wait=True)


def run_blocks(
    shape, block_chunks, brain_of, read, predict, out, tile, stride, n_workers
):
    """Predict every block that touches brain and write its votes into `out`.

    `shape` is that of the (c, z, y, x) inference array and `block_chunks`
    (from plan_chunks) lays out the blocks; `read(outer)` fetches one
    halo-grown block as a NumPy array.  `out` is a zarr array of the same
    shape whose chunks align with the blocks.  Blocks whose halo-grown extent
    has no brain are skipped without ever being read.  The rest are read and
    predicted from a pool of `n_workers` threads, so that while one block sits
    on a GPU the next ones are being fetched -- the reader (HDF5 for Imaris,
    zarr otherwise) rather than the GPU is the bottleneck for large inputs.
    `predict` serialises GPU access itself (see DevicePool).

    Reads are issued as they are needed rather than all up front, so after the
    first wave the pool settles into a pipeline: some threads reading, some
    on a GPU, rather than everyone reading at once and then everyone waiting.
    """
    halo = tile - stride
    plan = block_plan(shape[-3:], block_chunks[-3:], halo)
    todo = [item for item in plan if brain_of(item[1]).any()]
    print(
        f"{len(todo)} of {len(plan)} blocks touch brain; the rest are skipped unread",
        flush=True,
    )
    t0 = time.time()
    done = 0
    lock = threading.Lock()

    def work(item):
        nonlocal done
        core, outer, inner = item
        brain = brain_of(outer)
        vol = read(outer)  # (1, z, y, x): single input channel, asserted in main
        votes = predict(vol[0].astype(np.float32), brain)  # (n_out, z, y, x)
        out[(slice(None),) + core] = votes[(slice(None),) + inner]
        with lock:
            done += 1
            n = done
        elapsed = time.time() - t0
        print(
            f"block {n}/{len(todo)} done at {elapsed / 60:.1f} min "
            f"({elapsed / n:.0f} s/block, eta {(len(todo) - n) * elapsed / n / 60:.0f} min)",
            flush=True,
        )

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        list(pool.map(work, todo))


def main():
    tile = int(snakemake.params.tile)
    stride = int(snakemake.params.stride)
    batch_size = int(snakemake.params.batch_size)
    chunk = int(snakemake.params.chunk)
    n_gpus = int(snakemake.params.n_gpus)

    devices = resolve_devices(n_gpus)

    spec = dict(
        path=snakemake.input.spim,
        level=int(snakemake.params.start_level),
        target_scale=float(snakemake.params.iso_res),
        channel_labels=[snakemake.wildcards.stain],
        **snakemake.params.zarrnii_kwargs,
    )

    with get_dask_client("threads", snakemake.threads):
        znimg = load_near_isotropic(**spec)
        start_level = int(spec["level"])

        # Some acquisitions store a leading singleton time axis (t,c,z,y,x).
        # Inference runs on (c,z,y,x); the axis is restored on output so the
        # probseg keeps the same rank as its source store.
        has_t = znimg.data.ndim == 5 and znimg.data.shape[0] == 1
        data = open_inference_array(**spec)
        if data.ndim != 4 or data.shape[0] != 1:
            raise ValueError(
                f"expected a single-channel (c,z,y,x) or (1,c,z,y,x) image, "
                f"got shape {znimg.data.shape}"
            )
        block_chunks = plan_chunks(data.shape, chunk, tile - stride)

        brain_of = make_brain_sampler(snakemake.input.mask, data.shape)
        counter = TileCounter()

        pool = DevicePool(devices, snakemake.input.models)
        n_folds = len(pool.nets[devices[0]])
        n_classes = pool.nets[devices[0]][0].n_classes
        out_labels = ["plaque"] + (["falsepositive"] if n_classes >= 3 else [])
        res_um = {d: round(axis_scale_um(znimg, d), 2) for d in ("z", "y", "x")}
        print(
            f"loaded {n_folds} folds ({n_classes} classes) on each of "
            f"{[str(d) for d in devices]}; grid {data.shape} res {res_um} um "
            f"(from level {start_level}) chunk {chunk} halo {tile - stride} "
            f"tile {tile} stride {stride} batch {batch_size}",
            flush=True,
        )

        def predict(vol, brain_vol):
            with pool.acquire() as (device, nets):
                return predict_volume(
                    vol, nets, device, tile, stride, batch_size, brain_vol, counter
                )

        # Votes are staged in a temporary zarr, block by block, then read back
        # lazily for the pyramid writer below. uint8 and overwhelmingly zero, so
        # it compresses to a small fraction of the ~1 byte/voxel nominal size.
        tmp_dir = tempfile.TemporaryDirectory(
            dir=snakemake.resources.tmpdir, suffix="_lantern_votes"
        )
        tmp_store = os.path.join(tmp_dir.name, "votes.zarr")
        out = zarr.open_array(
            tmp_store,
            mode="w",
            shape=(len(out_labels),) + data.shape[1:],
            chunks=(1, chunk, chunk, chunk),
            dtype=np.uint8,
        )
        # Half the CPUs read (one process each, a couple of dask threads
        # inside), the rest of the budget is in-flight blocks waiting for or
        # sitting on a GPU.
        read = BlockReader(spec, n_procs=max(1, snakemake.threads // 2))
        try:
            run_blocks(
                data.shape,
                block_chunks,
                brain_of,
                read,
                predict,
                out,
                tile,
                stride,
                n_workers=snakemake.threads,
            )
        finally:
            read.close()
        votes = da.from_zarr(tmp_store).rechunk((len(out_labels),) + block_chunks[1:])

        # The fraction of folds voting each class: float32 in {0, .2, .4, .6, .8, 1}
        # for a 5-fold ensemble, one channel per class beyond background. Same
        # form as LANTERN's own whole-brain probmask.
        #
        # Deliberately a probability rather than the raw count or a 0-100 rescaling:
        # it is directly readable as model confidence in a viewer, and independent of
        # how many folds the ensemble happens to have. The 0/100 mask convention does
        # not apply -- that exists so fieldfrac can mean-pool a *mask* into a
        # percentage, and nothing but binarize_lantern_plaques reads this store.
        #
        # float32 costs ~4 bytes/voxel, so a level-1 whole brain is ~100 GB raw; it
        # compresses hard, being overwhelmingly zero.
        prob = votes.astype(np.float32) / float(n_folds)
        znimg_prob = with_channels(
            znimg, prob[np.newaxis] if has_t else prob, out_labels
        )

        with ProgressBar():
            znimg_prob.to_ome_zarr(
                snakemake.output.probseg,
                max_layer=5,
                match_scale_factors_from=snakemake.input.spim,
                **snakemake.config["zarrnii_out_kwargs"],
            )
        print(
            f"tiles run {counter.run}, skipped {counter.skipped} "
            f"(no brain in tile; blocks with no brain at all are not counted)",
            flush=True,
        )
        tmp_dir.cleanup()


if __name__ == "__main__":
    main()
