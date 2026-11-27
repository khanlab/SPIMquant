"""TEMPLATE: run a new deep-learning model through the tiled-inference engine.

Copy to ``spimquant/workflow/scripts/<yourmodel>_seg.py`` and fill in the
places marked ``TODO``.  Pick a script name that is NOT the name of a package
you import (a script called ``vesselfm.py`` shadows the ``vesselfm`` package
in tests and notebooks that put the scripts folder first on ``sys.path``).

The engine (``tiled_inference.py``) does the block planning, brain-mask
skipping, multi-process (Imaris-friendly) reading, overlapping-tile merging
and multi-GPU dispatch.  This file only says how to load the model and how to
turn a batch of raw tiles into per-tile outputs.  See
docs/howto/deep_learning_models.md for the full guide.

``tests/test_tiled_model_template.py`` runs this template with a tiny network,
so keep the class and function names if you change the template itself.
"""

import numpy as np
import torch
from dask.diagnostics import ProgressBar

from dask_setup import get_dask_client
from tiled_inference import ModelAdapter, TilingConfig, run_inference
from zarrnii import ZarrNii
from zarrnii_compat import with_channels


def zscore_tiles(tiles, eps=1e-8):
    """Per-tile z-score of a (B, z, y, x) batch, on the batch's device.

    TODO: replace with the normalisation the model was TRAINED with -- this is
    the most common source of silently bad output.  Do it per tile here (not
    over the whole image) unless the model was trained on whole-image
    statistics.  See vesselfm.scale_percentiles and
    lantern_plaques.normalize_tiles for percentile-based examples.
    """
    flat = tiles.reshape(tiles.shape[0], -1)
    mean = flat.mean(dim=1).view(-1, 1, 1, 1)
    std = flat.std(dim=1).view(-1, 1, 1, 1)
    return (tiles - mean) / std.clamp_min(eps)


class MyModel(ModelAdapter):
    """TODO: one line on what the model segments."""

    # One label per output channel; the stored image gets one channel each.
    out_labels = ["foreground"]  # TODO
    # dtype of predict_tiles' output.  Kept as-is under merge_mode="max";
    # the averaging modes ("gaussian", "average") always merge in float32.
    tile_dtype = np.float32

    def __init__(self, model_path, threshold=0.5):
        self.model_path = model_path
        self.threshold = float(threshold)

    def build_network(self):
        """The untrained network, matching the checkpoint's architecture.

        TODO: construct and return your torch.nn.Module, e.g. from the
        model's own package or a vendored definition (see lantern_model.py).
        """
        raise NotImplementedError("build the network here")

    def load(self, device):
        """Called once per device, from the main process, before any block.

        Return whatever predict_tiles needs (here the network); the engine
        hands it back unchanged.  Load from a file the workflow has already
        downloaded (an import_*_model rule) rather than from the network.
        """
        net = self.build_network()
        # TODO: adapt to the checkpoint layout (e.g. ckpt["state_dict"]).
        state = torch.load(self.model_path, map_location="cpu", weights_only=True)
        net.load_state_dict(state)
        return net.to(device).eval()

    def predict_tiles(self, net, tiles):
        """(B, t, t, t) float32 RAW tiles on the device -> (n_out, B, t, t, t).

        Runs on the device that `net` lives on; the block it came from holds
        that device exclusively, so no locking is needed.  Return a torch
        tensor (preferred: it stays on the device for the merge) or a NumPy
        array.
        """
        x = zscore_tiles(tiles)[:, None]  # (B, 1, t, t, t): one input channel
        with (
            torch.no_grad(),
            # bfloat16 is roughly 2x faster on A100s; only keep it if the
            # output matches fp32 on real data (compare on a crop first).
            torch.autocast(
                tiles.device.type,
                dtype=torch.bfloat16,
                enabled=tiles.device.type == "cuda",
            ),
        ):
            logits = net(x)  # TODO: (B, n_classes, t, t, t) for most U-Nets
        # Probabilities merge sensibly across overlapping tiles; logits and
        # hard labels do not.  TODO: pick the class channel(s) you output.
        prob = logits.float().sigmoid()[:, 0]  # (B, t, t, t)
        return prob[None]  # (n_out=1, B, t, t, t)

    def finalize(self, merged):
        """Merged (n_out, z, y, x) dask array -> what gets stored.

        Two conventions are used downstream:
        - a mask: uint8 0/100 (fieldfrac mean-pools it into a percentage);
        - a probability map: float32 in [0, 1], binarized by a later rule so
          the threshold can change without re-running inference (as LANTERN).
        """
        return (merged > self.threshold).astype(np.uint8) * np.uint8(100)

    def describe(self):
        return f"MyModel ({self.model_path}), threshold {self.threshold}"


def main():
    config = TilingConfig(
        tile=int(snakemake.params.tile),
        stride=int(snakemake.params.stride),
        batch_size=int(snakemake.params.batch_size),
        chunk=int(snakemake.params.chunk),
        merge_mode=str(snakemake.params.merge_mode),
        n_gpus=int(snakemake.params.n_gpus),
    )
    # What to read: the pyramid level in the output name, made
    # near-isotropic (as vesselFM).  For a model trained at a physical
    # resolution, give target_scale (um) instead and open `znimg` below with
    # zarrnii_compat.load_near_isotropic(**spec) (as LANTERN, whose output is
    # then named for that grid rather than for a level).
    spec = dict(
        path=snakemake.input.spim,
        level=int(snakemake.wildcards.level),
        target_scale=None,
        channel_labels=[snakemake.wildcards.stain],
        downsample_near_isotropic=True,
        **snakemake.params.zarrnii_kwargs,
    )
    adapter = MyModel(snakemake.input.model, snakemake.params.threshold)

    with get_dask_client("threads", snakemake.threads):
        # The image whose grid and metadata the output takes; opened the same
        # way as the engine opens `spec`.
        znimg = ZarrNii.from_file(
            **{k: v for k, v in spec.items() if k != "target_scale"}
        )
        has_t = znimg.data.ndim == 5 and znimg.data.shape[0] == 1

        with run_inference(
            adapter,
            spec,
            config,
            mask_path=snakemake.input.mask,
            threads=snakemake.threads,
            read_procs=snakemake.params.read_procs,
            tmpdir=snakemake.resources.tmpdir,
        ) as result:
            # `result` reads from a temporary store removed when the `with`
            # exits: write it out inside it.
            out = with_channels(
                znimg, result[np.newaxis] if has_t else result, adapter.out_labels
            )
            with ProgressBar():
                out.to_ome_zarr(
                    snakemake.output.mask,
                    max_layer=5,
                    match_scale_factors_from=snakemake.input.spim,
                    **snakemake.config["zarrnii_out_kwargs"],
                )


if __name__ == "__main__":
    main()
