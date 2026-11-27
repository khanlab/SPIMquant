"""Tiled deep-learning inference over large SPIM volumes, shared by every model.

A model plugs in through a small ``ModelAdapter`` (how to load it onto a device,
and how to turn a batch of raw tiles into per-tile outputs); everything else
lives here and is the same for every model:

- **Blocks.**  The volume is cut into ``chunk``-sized blocks (``plan_chunks``),
  each read with a halo of ``tile - stride`` voxels (``block_plan``) so that the
  tile overlap is honoured across block boundaries, not just within a block.
  Blocks are independent: each writes only its own core, so no cross-block
  merging is needed.

- **Brain-mask skipping.**  With a low-res brain mask, a block whose halo-grown
  extent has no brain is never read, and a tile with no brain is never run.
  Most of a light-sheet frame is air, so this is the main throughput lever.

- **Reading.**  Blocks are read in separate processes (``BlockReader``).  HDF5,
  which backs Imaris input, serialises every read in a process behind one
  global lock, so threads -- including dask's threaded scheduler -- cannot
  exceed one read stream; separate processes each get their own.

- **Tiling and merging.**  Within a block, ``tile``^3 tiles at ``stride`` are
  run in batches of ``batch_size`` and overlapping tiles are combined by
  ``merge_mode``: ``max``, ``average`` or ``gaussian`` (nnU-Net-style centre
  weighting, sigma = tile/8).

- **GPUs.**  The model is loaded once per device (``DevicePool``) and a block
  borrows a whole device for its duration, so several GPUs work in parallel
  and load-balance without a static split.  Reading, inference and writing
  run in a thread pool so that I/O overlaps the GPU.

``run_inference`` ties these together and yields the merged, adapter-finalised
output as a lazy dask array for the caller to write out.
"""

from __future__ import annotations

import contextlib
import itertools
import multiprocessing
import os
import queue
import tempfile
import threading
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
import torch

MERGE_MODES = ("average", "gaussian", "max")


# ----------------------------------------------------------------------
# model adapter


class ModelAdapter:
    """What the engine needs to know about one model (or ensemble).

    Subclasses set ``out_labels`` (one per output channel) and implement
    ``load`` and ``predict_tiles``; the rest have usable defaults.
    """

    #: one label per output channel; len(out_labels) is the channel count
    out_labels: list[str]
    #: dtype of predict_tiles' output, kept as-is under merge_mode="max"
    #: (averaging modes always accumulate in float32)
    tile_dtype = np.float32

    @property
    def n_out(self) -> int:
        return len(self.out_labels)

    def load(self, device: torch.device):
        """Return the model state for one device (called once per device)."""
        raise NotImplementedError

    def predict_tiles(self, model, tiles: torch.Tensor) -> torch.Tensor:
        """Per-tile outputs for one batch.

        ``tiles`` is (B, t, t, t) float32 RAW intensities, already on the
        model's device; preprocessing (normalisation, etc.) is the adapter's
        job.  Returns (n_out, B, t, t, t), as a tensor or array.
        """
        raise NotImplementedError

    def finalize(self, merged):
        """Map the merged (n_out, z, y, x) dask array to the stored output."""
        return merged.astype(np.float32)

    def describe(self) -> str:
        """One line for the run log."""
        return type(self).__name__


@dataclass
class TilingConfig:
    """Tiling, batching and device parameters, independent of the model.

    ``stride`` sets the overlap: ``tile - stride`` voxels between neighbouring
    tiles, which is also the halo read around each block.  ``from_overlap``
    builds one from a fractional overlap instead.
    """

    tile: int = 128
    stride: int = 64
    batch_size: int = 4
    chunk: int = 512
    merge_mode: str = "gaussian"
    n_gpus: int = 1

    def __post_init__(self):
        if self.merge_mode not in MERGE_MODES:
            raise ValueError(
                f"merge_mode must be one of {MERGE_MODES}, got {self.merge_mode!r}"
            )
        if not 0 < self.stride <= self.tile:
            raise ValueError(
                f"stride must be in (0, tile={self.tile}], got {self.stride}"
            )
        if self.chunk < self.tile:
            raise ValueError(f"chunk {self.chunk} is smaller than tile {self.tile}")

    @classmethod
    def from_overlap(cls, tile, overlap, **kwargs):
        """Config with stride from a fractional overlap in [0, 1)."""
        if not 0 <= overlap < 1:
            raise ValueError(f"overlap must be in [0, 1), got {overlap}")
        return cls(tile=tile, stride=max(1, round(tile * (1 - overlap))), **kwargs)

    @property
    def halo(self) -> int:
        return self.tile - self.stride


# ----------------------------------------------------------------------
# tiling and merging within one block


def tile_origins(extent, tile, stride):
    """Tile origins covering [0, extent), with the last pulled back to hit the edge."""
    if extent <= tile:
        return [0]
    origins = list(range(0, extent - tile + 1, stride))
    if origins[-1] + tile < extent:
        origins.append(extent - tile)
    return origins


def gaussian_importance_map(tile, sigma_scale=0.125, eps=1e-8):
    """Separable 3D Gaussian weight over one tile, 1 at its centre.

    sigma = tile * sigma_scale (tile/8, as in nnU-Net's sliding window).  The
    floor at `eps` keeps a tile's corners from weighing exactly zero, so a
    voxel covered only by tile corners (at the volume border) still gets a
    defined average.
    """
    coords = np.arange(tile, dtype=np.float32) - (tile - 1) / 2.0
    kernel = np.exp(-0.5 * (coords / (tile * sigma_scale)) ** 2).astype(np.float32)
    kernel /= kernel.max()
    weights = kernel[:, None, None] * kernel[None, :, None] * kernel[None, None, :]
    return np.maximum(weights, np.float32(eps))


def merge_tiles(acc, weight, batch_out, batch, tile, importance=None):
    """Fold one batch of per-tile outputs into a block's running result.

    `acc` is (n_out, z, y, x) and `batch_out` is (n_out, n_tiles, tile, tile,
    tile).  With `weight` None (merge_mode="max") the running per-voxel maximum
    is kept in `acc`.  Otherwise `acc` (float32) accumulates the weighted
    outputs and `weight` (z, y, x), shared by every channel, the weights --
    uniform when `importance` is None ("average"), else the Gaussian map --
    and predict_block divides one by the other at the end.
    """
    for j, (z, y, x) in enumerate(batch):
        window = (slice(z, z + tile), slice(y, y + tile), slice(x, x + tile))
        region = acc[(slice(None),) + window]
        if weight is None:
            np.maximum(region, batch_out[:, j], out=region)
        elif importance is None:
            region += batch_out[:, j]
            weight[window] += 1.0
        else:
            region += batch_out[:, j] * importance
            weight[window] += importance


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


def predict_block(
    vol,
    predict_tiles,
    *,
    n_out,
    device,
    tile,
    stride,
    batch_size,
    brain=None,
    counter=None,
    merge_mode="gaussian",
    tile_dtype=np.float32,
):
    """Tiled prediction over one 3D block, merged into (n_out, *vol.shape).

    `predict_tiles(tensor)` maps a (B, tile, tile, tile) float32 batch on
    `device` to (n_out, B, tile, tile, tile) outputs.  Overlapping tiles are
    combined per `merge_mode`: kept in `tile_dtype` for "max", float32
    (weighted) means for "average" and "gaussian".

    `brain`, if given, is a boolean array on the same grid as `vol`; tiles
    that contain no brain voxel are skipped and output zero.  Besides saving
    time, this keeps models whose per-tile normalisation would blow pure-air
    sensor noise up to unit variance from hallucinating far from the brain.

    The caller owns `device` exclusively for the duration of this call (see
    DevicePool), so no additional locking is needed here.
    """
    if merge_mode not in MERGE_MODES:
        raise ValueError(f"merge_mode must be one of {MERGE_MODES}, got {merge_mode!r}")
    orig_shape = vol.shape
    pad = [max(0, tile - n) for n in orig_shape]
    if any(pad):
        # Blocks at the volume border can be thinner than one tile; pad, then crop back.
        vol = np.pad(vol, [(0, p) for p in pad], mode="edge")
        if brain is not None:
            brain = np.pad(brain, [(0, p) for p in pad], mode="edge")

    if merge_mode == "max":
        acc = np.zeros((n_out,) + vol.shape, dtype=tile_dtype)
        weight = importance = None
    else:
        acc = np.zeros((n_out,) + vol.shape, dtype=np.float32)
        weight = np.zeros(vol.shape, dtype=np.float32)
        importance = gaussian_importance_map(tile) if merge_mode == "gaussian" else None
    crop = (slice(None),) + tuple(slice(0, n) for n in orig_shape)
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
        return acc[crop]

    for start in range(0, len(coords), batch_size):
        batch = coords[start : start + batch_size]
        arr = np.stack(
            [vol[z : z + tile, y : y + tile, x : x + tile] for z, y, x in batch]
        )
        out = predict_tiles(torch.from_numpy(arr).to(device))
        if isinstance(out, torch.Tensor):
            out = out.cpu().numpy()
        merge_tiles(acc, weight, out, batch, tile, importance)

    if weight is not None:
        # Voxels that no tile covered (every covering tile was skipped for having
        # no brain) keep weight 0 and stay at 0.
        np.divide(acc, weight, out=acc, where=weight > 0)
    return acc[crop]


# ----------------------------------------------------------------------
# devices


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
    """Lends out one device at a time, with the model resident on each.

    `load(device)` is called once per device.  Handing out devices from a
    queue rather than striping blocks by index balances the load for free:
    blocks vary enormously in cost because most of a frame is background, so a
    static assignment would leave one card idle.
    """

    def __init__(self, devices, load):
        self._free = queue.Queue()
        self.models = {}
        for device in devices:
            self.models[device] = load(device)
            self._free.put(device)

    @property
    def devices(self):
        return list(self.models)

    @contextlib.contextmanager
    def acquire(self):
        device = self._free.get()
        try:
            yield device, self.models[device]
        finally:
            self._free.put(device)


# ----------------------------------------------------------------------
# blocks


def make_brain_sampler(mask_path, image_shape):
    """A function mapping (z, y, x) slices on the inference grid to a brain mask.

    With `mask_path` None every voxel counts as brain.  Otherwise the mask
    (a NIfTI at a much coarser level) is held in memory and nearest-neighbour
    sampled into whatever region is asked for, from the region's array
    location.  It is first dilated by one low-res voxel (tens of microns on
    the inference grid) so that a tile at the pial surface is never skipped on
    account of a slightly tight mask.
    """
    if mask_path is None:
        return lambda slices: np.ones(
            tuple(sl.stop - sl.start for sl in slices), dtype=bool
        )

    from scipy.ndimage import binary_dilation
    from zarrnii import ZarrNii

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
    halo (z=1030 at block 512 leaves a 6-voxel tail), and a block that thin
    cannot be bridged by the halo of its neighbour.  A short tail is absorbed
    into its predecessor instead.  Every chunk start stays a multiple of
    `block`, so the per-block tile grid stays aligned with the global one.
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
    lying wholly inside outer -- which is what makes the tile overlap correct
    across block boundaries and not just within them.  An axis that is a
    single chunk gets no halo: there is no block boundary to bridge.
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
    from zarrnii_compat import open_inference_array

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
    """Predict every block that touches brain and write its output into `out`.

    `shape` is that of the (c, z, y, x) inference array and `block_chunks`
    (from plan_chunks) lays out the blocks; `read(outer)` fetches one
    halo-grown block as a NumPy array.  `out` is a zarr array of the same
    spatial shape whose chunks align with the blocks.  Blocks whose halo-grown
    extent has no brain are skipped without ever being read.  The rest are
    read and predicted from a pool of `n_workers` threads, so that while one
    block sits on a GPU the next ones are being fetched -- the reader rather
    than the GPU is usually the bottleneck for large inputs.  `predict`
    serialises GPU access itself (see DevicePool).

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
        vol = read(outer)  # (1, z, y, x): single input channel, asserted by caller
        result = predict(vol[0].astype(np.float32), brain)  # (n_out, z, y, x)
        out[(slice(None),) + core] = result[(slice(None),) + inner]
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


# ----------------------------------------------------------------------
# driver


@contextlib.contextmanager
def run_inference(adapter, spec, config, *, mask_path=None, threads=8, tmpdir=None):
    """Run `adapter` over the image described by `spec`; yield the output.

    `spec` holds the arguments to ``zarrnii_compat.open_inference_array``
    (path, level, target_scale, channel_labels and any from_file kwargs), and
    must select a single channel.  Yields a lazy (n_out, z, y, x) dask array,
    already passed through ``adapter.finalize`` and chunked like the blocks;
    it reads from a temporary zarr that is removed when the context exits, so
    write it out inside the ``with``.

    `threads` is the CPU budget: half go to reader processes, all of it to
    in-flight blocks (reading, waiting for a device, or on one).
    """
    import dask.array as da
    import zarr
    from zarrnii_compat import open_inference_array

    data = open_inference_array(**spec)
    if data.ndim != 4 or data.shape[0] != 1:
        raise ValueError(
            f"expected a single-channel (c,z,y,x) image, got shape {data.shape}"
        )
    block_chunks = plan_chunks(data.shape, config.chunk, config.halo)
    brain_of = make_brain_sampler(mask_path, data.shape)
    counter = TileCounter()

    devices = resolve_devices(config.n_gpus)
    pool = DevicePool(devices, adapter.load)
    print(
        f"{adapter.describe()} on {[str(d) for d in devices]}; grid {data.shape} "
        f"chunk {config.chunk} halo {config.halo} tile {config.tile} "
        f"stride {config.stride} batch {config.batch_size} merge {config.merge_mode}",
        flush=True,
    )

    def predict(vol, brain_vol):
        with pool.acquire() as (device, model):
            return predict_block(
                vol,
                lambda tiles: adapter.predict_tiles(model, tiles),
                n_out=adapter.n_out,
                device=device,
                tile=config.tile,
                stride=config.stride,
                batch_size=config.batch_size,
                brain=brain_vol,
                counter=counter,
                merge_mode=config.merge_mode,
                tile_dtype=adapter.tile_dtype,
            )

    # Outputs are staged in a temporary zarr, block by block, then read back
    # lazily for the caller's writer.  Overwhelmingly zero outside the brain,
    # so it compresses to a small fraction of its nominal size.
    with tempfile.TemporaryDirectory(dir=tmpdir, suffix="_tiled_inference") as tmp:
        store = os.path.join(tmp, "merged.zarr")
        out = zarr.open_array(
            store,
            mode="w",
            shape=(adapter.n_out,) + data.shape[1:],
            chunks=(1,) + (config.chunk,) * 3,
            dtype=adapter.tile_dtype if config.merge_mode == "max" else np.float32,
        )
        read = BlockReader(spec, n_procs=max(1, threads // 2))
        try:
            run_blocks(
                data.shape,
                block_chunks,
                brain_of,
                read,
                predict,
                out,
                config.tile,
                config.stride,
                n_workers=threads,
            )
        finally:
            read.close()
        print(
            f"tiles run {counter.run}, skipped {counter.skipped} "
            f"(no brain in tile; blocks with no brain at all are not counted)",
            flush=True,
        )
        merged = da.from_zarr(store).rechunk((adapter.n_out,) + block_chunks[1:])
        yield adapter.finalize(merged)
