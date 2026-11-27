"""5-fold LANTERN ensemble Abeta plaque segmentation.

Runs a published LANTERN ensemble (``apooladi/lantern-<name>``, chosen with
``--lantern_config``; ``abeta`` by default) over the Abeta channel and
writes the per-fold agreement map -- the fraction of folds voting foreground,
averaged over the overlapping tiles that saw each voxel (``{0, .2, .4, .6, .8,
1}`` for five folds only under ``merge_mode="max"``) -- as an OME-Zarr.
Binarization is deliberately left to a downstream rule so the vote threshold
can be changed without re-running inference.

Every fold is a decoder fine-tuned on the same frozen self-supervised encoder,
so the encoder is run once per batch and its skip features are fed to each
fold's decoder.  The encoder is ~2/3 of the FLOPs of a full forward pass, so
this roughly halves inference time against five independent passes, and it is
exact: the encoder is deterministic in eval mode and the decoders only read the
skips.  ``load_folds`` checks the encoders really are identical rather than
trusting the model card.  An ensemble with a bottleneck conditioner (below)
keeps that saving: the conditioner is trained with the decoder, so it is
per-fold and runs inside ``ResEncLUNet.decode``, after the one shared encoder
pass.

The block planning, brain-mask skipping, multi-process reading, tiling,
overlap merging and multi-GPU dispatch are model-independent and live in
``tiled_inference``; this script supplies only the ``LanternEnsemble``
adapter (loading the folds, the per-tile normalisation and the fold votes)
and the output conventions.  Notes specific to LANTERN:

- The model reads RAW intensities.  It was never trained on N4-corrected or
  otherwise rescaled data, so this path takes ``input.spim`` directly rather than
  the ``desc-corrected*`` store the GMM/Otsu rules consume.

- The model is scale-sensitive: it was trained at ~4 um isotropic.  The input
  is loaded via ``load_near_isotropic`` rather than at a fixed pyramid level,
  because a fixed level lands on the right grid only for pyramids that
  downsample x/y alone (ki3); a pyramid that also halves z per level (mapt)
  puts z at 8 um by level 1, so the loader walks to a finer level and
  downsamples per axis to get near ``plaque_iso_res``.

Inference contract (from the model cards, must match training exactly):
  tiles of 128^3 at stride 64, per-tile 0.5/99.5 percentile clip then z-score,
  softmax over the classes, plaque where class-1 probability >= 0.5 per fold.
  3-class ensembles (classes background / plaque / "annotated false positive")
  additionally emit, as a second output channel, the number of folds whose
  argmax is class 2 -- the model card's rule is "class 2 where any fold's
  argmax is class 2", which is that channel at >= 1 vote, and
  binarize_lantern_plaques applies it.

  Ensembles carrying a bottleneck conditioner (``apooladi/lantern-abeta``, the
  default) additionally require the mean and standard deviation of each RAW
  tile, i.e. the two numbers the clip-and-z-score above destroys; they are fed
  to the bottleneck as extra features (see lantern_model.BottleneckConditioner)
  and the forward pass refuses to run without them.  ``raw_tile_stats``
  measures them on the tile as read, before normalisation, which is how the
  standardisation constants in the checkpoints were measured.  Whether an
  ensemble wants them is read from the checkpoints, not configured, so the
  older unconditioned ensembles (``ki3-abeta``, ``mapt_a/b/c-abeta``) keep
  running unchanged.

Overlapping tiles (up to 8 per voxel at stride = tile/2) are combined by
``merge_mode``, the one deliberate departure from the model cards:

- ``max`` (the published LANTERN inference) keeps each voxel's highest
  per-tile vote count.  It recovers a plaque clipped at one tile's edge, but a
  tile that saw the voxel at its own edge, with little context (e.g. at the
  brain surface), can vote 5/5 on its own and win -- the tiling artifacts.
- ``average`` takes the mean vote count over the tiles that ran, so that tile
  is outvoted by the rest.
- ``gaussian`` (default) weights each tile by a Gaussian centred on it
  (sigma = tile/8, as in nnU-Net) before averaging, so tile edges count least
  -- the same idea as ``average``, but it discounts a tile's edge rather than
  trusting it as much as its centre.

The averaged modes give fractional counts; binarize_lantern_plaques cuts them
at the midpoint below the requested vote threshold (2.5 for "3 of 5").
"""

import numpy as np
import torch
from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from lantern_model import ResEncLUNet
from tiled_inference import (
    MERGE_MODES,
    ModelAdapter,
    TilingConfig,
    gaussian_importance_map,
    merge_tiles,
    predict_block,
    run_inference,
)
from zarrnii_compat import axis_scale_um, load_near_isotropic, with_channels

# Former name, kept for the tests and anything else that imports it from here
# (as are MERGE_MODES and gaussian_importance_map, imported above).
merge_tile_votes = merge_tiles


def raw_tile_stats(x):
    """Per-tile mean and std of the RAW intensities, as [B, 1] conditioner inputs.

    Must be called on the batch BEFORE normalize_tiles: these are exactly the
    two numbers that transform removes, which is why the conditioner exists.
    `correction=0` is numpy's `.std()` convention, the one the checkpoints'
    standardisation constants were measured with (over the raw 128^3 training
    patches); at 2M voxels per tile the difference from the unbiased estimator
    is ~1e-7 relative, but matching it costs nothing.
    """
    flat = x.reshape(x.shape[0], -1).float()
    return (
        flat.mean(dim=1, keepdim=True),
        flat.std(dim=1, correction=0, keepdim=True),
    )


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


def vote_tiles(nets, tensor):
    """Per-tile fold vote counts for one batch of RAW tiles, as uint8.

    Returns (n_out, B, t, t, t): channel 0 counts the folds that call plaque
    (class-1 probability >= 0.5); for a 3-class ensemble channel 1 counts the
    folds whose argmax is class 2.

    Folds are summed WITHIN a tile here and tiles are combined afterwards by
    the engine's merge_mode.  Under "max" those two steps do not commute --
    max-then-sum would score a voxel higher when different folds find it from
    different tiles -- but this order is the one LANTERN's own inference uses
    (scripts/infer_wholebrain.py), and plaque_vote_threshold was calibrated
    under it.
    """
    device = tensor.device
    n_out = 2 if nets[0].n_classes >= 3 else 1
    tile = tensor.shape[1:]
    # Before normalize_tiles: the conditioner wants the raw tile's moments.
    raw_mean, raw_std = raw_tile_stats(tensor)
    tensor = normalize_tiles(tensor)[:, None]
    out = torch.zeros((n_out, tensor.shape[0], *tile), dtype=torch.uint8, device=device)
    with (
        torch.no_grad(),
        torch.autocast(
            device.type, dtype=torch.bfloat16, enabled=device.type == "cuda"
        ),
    ):
        # One encoder pass serves every fold: load_folds has verified the
        # encoders are identical and pointed all folds at the same module.
        # Each fold's decode() then applies its OWN bottleneck conditioner,
        # if it has one, to a copy of the deepest skip.
        skips = nets[0].encoder(tensor)
        for net in nets:
            logits = net.decode(skips, raw_mean, raw_std)
            prob = torch.softmax(logits.float(), dim=1)
            out[0] += (prob[:, 1] >= 0.5).to(torch.uint8)
            if n_out == 2:
                out[1] += (prob.argmax(dim=1) == 2).to(torch.uint8)
    return out


def predict_volume(
    vol,
    nets,
    device,
    tile,
    stride,
    batch_size,
    brain=None,
    counter=None,
    merge_mode="gaussian",
):
    """Tiled fold-vote maps for one 3D block, as counts in [0, len(nets)].

    Returns (n_out, *vol.shape): uint8 counts for "max", float32 (weighted)
    mean counts for "average" and "gaussian".  See
    tiled_inference.predict_block for the tiling and brain skipping.
    """
    return predict_block(
        vol,
        lambda tiles: vote_tiles(nets, tiles),
        n_out=2 if nets[0].n_classes >= 3 else 1,
        device=device,
        tile=tile,
        stride=stride,
        batch_size=batch_size,
        brain=brain,
        counter=counter,
        merge_mode=merge_mode,
        tile_dtype=np.uint8,
    )


def share_encoder(nets, model_paths):
    """Point every fold at fold 0's encoder module, after checking they match.

    LANTERN fine-tunes only the decoder (and, where present, the bottleneck
    conditioner), on an encoder frozen from self-supervised pretraining, so
    every fold of an ensemble carries the same encoder weights.  Only the
    encoder is shared: the conditioner's weights and standardisation constants
    are fitted per fold on that fold's own training split, so each fold keeps
    its own.  That is what lets ``predict_volume`` run the encoder once
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
    tiled_inference.DevicePool then guarantees only one block uses a given
    device at a time, which bounds peak VRAM without needing a lock inside the forward pass.

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
        cond_dim = checkpoint_cond_dim(ckpt, path)
        net = ResEncLUNet(
            num_classes=n_classes,
            num_input_channels=1,
            conditioner=cond_dim is not None,
            cond_dim=cond_dim or 32,
        )
        net.load_state_dict(ckpt["state_dict"])
        net.n_classes = n_classes
        nets.append(net.eval())
    if len({net.n_classes for net in nets}) != 1:
        raise ValueError(
            f"folds disagree on the number of classes: "
            f"{[net.n_classes for net in nets]}"
        )
    # A conditioned fold needs the raw statistics and an unconditioned one
    # ignores them, so a mixed ensemble would silently weight two different
    # models together rather than fail.
    if len({net.cond is not None for net in nets}) != 1:
        raise ValueError(
            "folds disagree on whether they have a bottleneck conditioner: "
            + ", ".join(
                f"{path}={'yes' if net.cond is not None else 'no'}"
                for net, path in zip(nets, model_paths)
            )
        )
    share_encoder(nets, list(model_paths))
    return [net.to(device) for net in nets]


def checkpoint_cond_dim(ckpt, path):
    """`cond_dim` of a checkpoint's bottleneck conditioner, or None if it has none.

    Read from the state_dict rather than from the metadata, because the
    state_dict is what has to load: `cond.proj` is Linear(n_cond, cond_dim), so
    its weight gives both numbers.  The `bottleneck_cond` metadata block is
    cross-checked against it so a checkpoint whose card and weights disagree is
    an error rather than a silently unconditioned run.

    n_cond must be 2 (log raw mean, log raw std of one stain).  Upstream also
    trains 4-slot multi-stain conditioners, whose 12 or 16 inputs this
    single-channel pipeline cannot supply.
    """
    weight = ckpt["state_dict"].get("cond.proj.weight")
    declared = ckpt.get("bottleneck_cond") is not None
    if (weight is not None) != declared:
        raise ValueError(
            f"{path}: bottleneck_cond metadata says "
            f"{'conditioned' if declared else 'unconditioned'} but the "
            f"state_dict {'has' if weight is not None else 'has no'} cond.proj"
        )
    if weight is None:
        return None
    cond_dim, n_cond = int(weight.shape[0]), int(weight.shape[1])
    if n_cond != 2:
        raise ValueError(
            f"{path}: conditioner takes {n_cond} inputs; this pipeline feeds a "
            "single stain and can only supply 2 (log raw mean, log raw std)"
        )
    return cond_dim


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


class LanternEnsemble(ModelAdapter):
    """A LANTERN ensemble of any number of folds (one is fine) as an adapter.

    The stored output is the fraction of folds voting each class, so it reads
    the same whatever the fold count.
    """

    tile_dtype = np.uint8  # fold counts

    def __init__(self, model_paths):
        self.model_paths = list(model_paths)
        self.n_folds = len(self.model_paths)

    def load(self, device):
        # The layout is read from the checkpoints, so it is only known once the
        # first device has loaded them; every device loads the same files.
        nets = load_folds(self.model_paths, device)
        self.n_classes = nets[0].n_classes
        self.conditioned = nets[0].cond is not None
        self.out_labels = ["plaque"] + (
            ["falsepositive"] if self.n_classes >= 3 else []
        )
        return nets

    def predict_tiles(self, nets, tiles):
        return vote_tiles(nets, tiles)

    def finalize(self, merged):
        # The fraction of folds voting each class: float32 in [0, 1], and in
        # {0, .2, .4, .6, .8, 1} for a 5-fold ensemble under max merging (the
        # same form as LANTERN's own whole-brain probmask).
        #
        # Deliberately a probability rather than the raw count or a 0-100
        # rescaling: it is directly readable as model confidence in a viewer,
        # and independent of how many folds the ensemble happens to have. The
        # 0/100 mask convention does not apply -- that exists so fieldfrac can
        # mean-pool a *mask* into a percentage, and nothing but
        # binarize_lantern_plaques reads this store.
        return merged.astype(np.float32) / float(self.n_folds)

    def describe(self):
        return (
            f"LANTERN {self.n_folds} fold(s), {self.n_classes} classes, "
            f"{'with' if self.conditioned else 'no'} bottleneck conditioner"
        )


def main():
    config = TilingConfig(
        tile=int(snakemake.params.tile),
        stride=int(snakemake.params.stride),
        batch_size=int(snakemake.params.batch_size),
        chunk=int(snakemake.params.chunk),
        merge_mode=str(snakemake.params.merge_mode),
        n_gpus=int(snakemake.params.n_gpus),
    )
    spec = dict(
        path=snakemake.input.spim,
        level=int(snakemake.params.start_level),
        target_scale=float(snakemake.params.iso_res),
        channel_labels=[snakemake.wildcards.stain],
        **snakemake.params.zarrnii_kwargs,
    )
    adapter = LanternEnsemble(snakemake.input.models)

    with get_dask_client("threads", snakemake.threads):
        znimg = load_near_isotropic(**spec)
        res_um = {d: round(axis_scale_um(znimg, d), 2) for d in ("z", "y", "x")}
        print(f"res {res_um} um (from level {spec['level']})", flush=True)
        # Some acquisitions store a leading singleton time axis (t,c,z,y,x).
        # Inference runs on (c,z,y,x); the axis is restored on output so the
        # probseg keeps the same rank as its source store.
        has_t = znimg.data.ndim == 5 and znimg.data.shape[0] == 1

        with run_inference(
            adapter,
            spec,
            config,
            mask_path=snakemake.input.mask,
            threads=snakemake.threads,
            read_procs=snakemake.params.read_procs,
            tmpdir=snakemake.resources.tmpdir,
        ) as prob:
            # float32 costs ~4 bytes/voxel, so a level-1 whole brain is ~100 GB
            # raw; it compresses hard, being overwhelmingly zero.
            znimg_prob = with_channels(
                znimg, prob[np.newaxis] if has_t else prob, adapter.out_labels
            )
            with ProgressBar():
                znimg_prob.to_ome_zarr(
                    snakemake.output.probseg,
                    max_layer=5,
                    **snakemake.config["zarrnii_out_kwargs"],
                )


if __name__ == "__main__":
    main()
