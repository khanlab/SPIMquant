"""vesselFM vessel segmentation, run through the shared tiled-inference engine.

The network, checkpoint and preprocessing are vesselFM's own (from the
``vesselfm`` package's bundled ``inference_zarrnii`` config); the tiling,
overlap merging, brain-mask skipping, multi-process Imaris reads and multi-GPU
dispatch come from ``tiled_inference`` (see that module).

Differences from the former ``ZarrNii.segment(VesselFMPlugin)`` path, which
ran each disjoint 128^3 dask block as a single patch:

- Tiles now overlap (``vessel_stride``, half a tile by default, as in
  vesselFM's own config) and are merged per ``vessel_merge_mode``, including
  across block boundaries, so there are no seams at block edges.  What is
  merged is the per-tile sigmoid probability; the threshold is applied after
  merging.
- Preprocessing (1st/99th-percentile rescale to [0, 1], clipped) is still per
  128^3 tile, as before, but runs on the GPU over a whole batch.
- Tiles and blocks without brain in the low-res brain mask are skipped and
  stay 0.
"""

import numpy as np
import torch
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii

from dask_setup import get_dask_client
from tiled_inference import ModelAdapter, TilingConfig, run_inference
from zarrnii_compat import with_channels


def scale_percentiles(tiles, lower=0.01, upper=0.99, eps=1e-8):
    """Per-tile percentile rescale to [0, 1], clipped.

    The torch equivalent of vesselFM's MONAI ``ScaleIntensityRangePercentiles
    (lower=1, upper=99, b_min=0, b_max=1, clip=True)``, applied to each tile of
    a (B, z, y, x) batch separately.  Both use linear interpolation between
    order statistics.
    """
    flat = tiles.reshape(tiles.shape[0], -1).float()
    q = torch.tensor([lower, upper], device=tiles.device, dtype=flat.dtype)
    lo, hi = torch.quantile(flat, q, dim=1)  # each (B,)
    shape = (-1,) + (1,) * (tiles.ndim - 1)
    lo, hi = lo.view(shape), hi.view(shape)
    return ((tiles - lo) / (hi - lo).clamp_min(eps)).clamp(0.0, 1.0)


class VesselFM(ModelAdapter):
    """vesselFM as an adapter: per-tile sigmoid vessel probability."""

    out_labels = ["vessels"]
    tile_dtype = np.float32

    def __init__(self, model_path, threshold=0.5):
        self.model_path = model_path
        self.threshold = float(threshold)

    def load(self, device):
        from vesselfm.zarrnii_plugin.vesselfm import create_config, load_model

        model = load_model(create_config(), device, model_path=self.model_path)
        return model.to(device).eval()

    def predict_tiles(self, model, tiles):
        x = scale_percentiles(tiles)[:, None]  # (B, 1, t, t, t)
        with torch.no_grad():
            logits = model(x)
        return logits[:, 0].float().sigmoid()[None]  # (1, B, t, t, t)

    def finalize(self, merged):
        # 0/100 mask, the convention fieldfrac mean-pools into a percentage.
        return (merged > self.threshold).astype(np.uint8) * np.uint8(100)

    def describe(self):
        return f"vesselFM ({self.model_path}), threshold {self.threshold}"


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
        level=int(snakemake.wildcards.level),
        target_scale=None,
        channel_labels=[snakemake.wildcards.stain],
        downsample_near_isotropic=True,
        **snakemake.params.zarrnii_kwargs,
    )
    adapter = VesselFM(snakemake.input.model_path, snakemake.params.threshold)

    with get_dask_client("threads", snakemake.threads):
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
        ) as mask:
            znimg_mask = with_channels(
                znimg, mask[np.newaxis] if has_t else mask, [snakemake.wildcards.stain]
            )
            with ProgressBar():
                znimg_mask.to_ome_zarr(
                    snakemake.output.mask,
                    max_layer=5,
                    match_scale_factors_from=snakemake.input.spim,
                    **snakemake.config["zarrnii_out_kwargs"],
                )


if __name__ == "__main__":
    main()
