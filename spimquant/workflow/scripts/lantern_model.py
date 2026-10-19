"""ResEncL U-Net architecture for the LANTERN Abeta plaque segmentation ensemble.

Vendored from the LANTERN repository (``lantern/encoders/resenc_l.py``,
``lantern/models/unet_seg.py`` and ``lantern/models/slot_stem.py``) so that
SPIMquant can load the published checkpoints without depending on that
repository.  The definitions must stay byte-compatible with the checkpoints on
HuggingFace (``apooladi/lantern-abeta`` and its cohort-specific siblings) --
changing any architectural constant here will make ``load_state_dict`` fail.

The encoder is identical to nnssl's ``architecture_registry.get_res_enc_l()``;
the decoder is nnU-Net's ``UNetDecoder`` over that encoder's skips.

Only the single-channel path is vendored.  Upstream's ``ResEncLUNet`` can also
put a 4-slot multi-stain stem in front of the encoder and run several output
heads; no published Abeta ensemble uses either, and SPIMquant feeds the model
one stain, so ``build_resenc_l`` stands in for upstream's ``build_backbone``
and there is one head.
"""

from __future__ import annotations

import torch
import torch.nn as nn
from dynamic_network_architectures.building_blocks.residual_encoders import (
    ResidualEncoder,
)
from dynamic_network_architectures.building_blocks.unet_decoder import UNetDecoder

N_COORDS = 3  # upstream's (AP, DV, |ML|) atlas coordinate, unused here


def build_resenc_l(num_input_channels: int = 1) -> ResidualEncoder:
    """Build the ResEncL encoder.

    Architecture: 6 stages, features [32, 64, 128, 256, 320, 320],
    strides [1,1,1] then [2,2,2] x 5, blocks per stage [1, 3, 4, 6, 6, 6],
    InstanceNorm3d (eps=1e-5, affine=True), LeakyReLU(inplace=True).
    """
    n_stages = 6
    encoder = ResidualEncoder(
        input_channels=num_input_channels,
        n_stages=n_stages,
        features_per_stage=[32, 64, 128, 256, 320, 320],
        conv_op=nn.Conv3d,
        kernel_sizes=[[3, 3, 3]] * n_stages,
        strides=[[1, 1, 1]] + [[2, 2, 2]] * 5,
        n_blocks_per_stage=[1, 3, 4, 6, 6, 6],
        conv_bias=True,
        norm_op=nn.InstanceNorm3d,
        norm_op_kwargs={"eps": 1e-5, "affine": True},
        nonlin=nn.LeakyReLU,
        nonlin_kwargs={"inplace": True},
        return_skips=True,
        disable_default_stem=False,
        stem_channels=None,
    )
    return encoder


class BottleneckConditioner(nn.Module):
    """Feed the bottleneck the absolute brightness that per-tile normalisation removes.

    LANTERN's input transform clips and z-scores every tile, which throws away
    the tile's absolute mean and spread -- information a plaque segmenter wants,
    since the same z-scored texture means different things in a dim and a bright
    field.  The conditioner gives those two numbers back: standardised
    ``log1p(raw_mean)`` and ``log1p(raw_std)`` of the RAW tile go through
    ``Linear(2 -> cond_dim)``, broadcast to every voxel of the bottleneck feature
    map, and are concatenated and mixed back to ``channels`` by a 1x1x1 conv.

    The standardisation constants ride in the checkpoint as buffers, measured on
    that fold's own training split, so each fold of an ensemble carries its own
    (unlike the frozen encoder, which they share).  ``mix`` is identity-initialised
    -- identity on the features, zero on the conditioning -- so switching the
    conditioner on cannot change a warm-started network at step 0; a loaded
    checkpoint overwrites it either way.

    Upstream also supports per-slot statistics for its 4-slot multi-stain input,
    an optional Allen atlas coordinate, and dropout on the two statistics during
    training.  SPIMquant runs one stain, passes no coordinate and only ever runs
    in eval mode, so this copy hardcodes ``n_slots=1`` with no present flag (it
    would be a constant 1) and no coordinate.  The coordinate buffers are still
    registered because upstream registers them unconditionally and they are
    therefore in the published state_dicts.
    """

    def __init__(self, channels: int = 320, cond_dim: int = 32) -> None:
        super().__init__()
        self.n_cond = 2  # log raw mean, log raw std, of the one stain
        self.proj = nn.Linear(self.n_cond, cond_dim)
        self.mix = nn.Conv3d(channels + cond_dim, channels, kernel_size=1, bias=True)
        with torch.no_grad():
            self.mix.weight.zero_()
            for c in range(channels):
                self.mix.weight[c, c, 0, 0, 0] = 1.0
            self.mix.bias.zero_()
        # Standardisation of the two log statistics; identity until a checkpoint loads.
        for name in ("log_mean_mu", "log_mean_sd", "log_std_mu", "log_std_sd"):
            self.register_buffer(
                name, torch.zeros(1) if name.endswith("mu") else torch.ones(1)
            )
        # Upstream's atlas-coordinate frame: never read here, present in the state_dict.
        self.register_buffer("coord_axis", torch.arange(N_COORDS, dtype=torch.long))
        self.register_buffer("coord_center", torch.zeros(N_COORDS))
        self.register_buffer("coord_half", torch.ones(N_COORDS))
        self.register_buffer("coord_abs", torch.tensor([0.0, 0.0, 1.0]))
        self.register_buffer("stats_loaded", torch.zeros(1))

    def cond_vector(
        self, raw_mean: torch.Tensor, raw_std: torch.Tensor
    ) -> torch.Tensor:
        """[B, 2]: the standardised log raw mean and log raw std.

        `raw_mean` and `raw_std` are [B, 1] statistics of the raw, unnormalised
        tile.  `clamp_min(0)` guards log1p against a negative mean, which a
        background-subtracted store could produce.
        """
        lm = (torch.log1p(raw_mean.clamp_min(0)) - self.log_mean_mu) / self.log_mean_sd
        ls = (torch.log1p(raw_std.clamp_min(0)) - self.log_std_mu) / self.log_std_sd
        return torch.cat([lm, ls], dim=1)

    def forward(
        self, feat: torch.Tensor, raw_mean: torch.Tensor, raw_std: torch.Tensor
    ) -> torch.Tensor:
        """`feat` is the bottleneck [B, C, d, h, w]; returns the same shape."""
        v = self.cond_vector(raw_mean, raw_std)
        c = self.proj(v.to(feat.dtype))[:, :, None, None, None].expand(
            -1, -1, *feat.shape[2:]
        )
        return self.mix(torch.cat([feat, c], dim=1))


class ResEncLUNet(nn.Module):
    """ResEncL encoder + UNetDecoder, 2- or 3-class (background / plaque / FP) logits.

    With `conditioner=True` a BottleneckConditioner sits on the deepest skip, and
    the raw tile statistics become required inputs to the forward pass.
    """

    def __init__(
        self,
        num_classes: int = 2,
        num_input_channels: int = 1,
        n_conv_per_stage_decoder: int = 1,
        deep_supervision: bool = False,
        conditioner: bool = False,
        cond_dim: int = 32,
    ) -> None:
        super().__init__()
        self.encoder = build_resenc_l(num_input_channels)
        n_dec_stages = len(self.encoder.output_channels) - 1
        self.decoder = UNetDecoder(
            self.encoder,
            num_classes,
            [n_conv_per_stage_decoder] * n_dec_stages,
            deep_supervision=deep_supervision,
        )
        self.cond = (
            BottleneckConditioner(
                channels=self.encoder.output_channels[-1], cond_dim=cond_dim
            )
            if conditioner
            else None
        )

    def decode(
        self,
        skips,
        raw_mean: torch.Tensor | None = None,
        raw_std: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Decode encoder `skips`, conditioning the bottleneck if this fold has a conditioner.

        Split out of `forward` so that an ensemble sharing one frozen encoder can
        run it once per batch and still condition per fold -- the conditioner is
        fine-tuned alongside the decoder, so it is NOT shared (see
        lantern_plaques.share_encoder).  `skips` is copied rather than written
        through, because the caller hands the same list to every fold.
        """
        if self.cond is not None:
            if raw_mean is None or raw_std is None:
                raise ValueError(
                    "the bottleneck conditioner needs raw_mean and raw_std"
                )
            skips = list(skips)
            skips[-1] = self.cond(skips[-1], raw_mean, raw_std)
        return self.decoder(skips)

    def forward(
        self,
        x: torch.Tensor,
        raw_mean: torch.Tensor | None = None,
        raw_std: torch.Tensor | None = None,
    ) -> torch.Tensor:
        return self.decode(self.encoder(x), raw_mean, raw_std)
