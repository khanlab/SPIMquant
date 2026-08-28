"""ResEncL U-Net architecture for the LANTERN Abeta plaque segmentation ensemble.

Vendored from the LANTERN repository (``lantern/encoders/resenc_l.py`` and
``lantern/models/unet_seg.py``) so that SPIMquant can load the published
checkpoints without depending on that repository.  The definitions must stay
byte-compatible with the checkpoints on HuggingFace
(``apooladi/lantern-ki3-abeta``) -- changing any architectural constant here
will make ``load_state_dict`` fail.

The encoder is identical to nnssl's ``architecture_registry.get_res_enc_l()``;
the decoder is nnU-Net's ``UNetDecoder`` over that encoder's skips.
"""

from __future__ import annotations

import torch
import torch.nn as nn
from dynamic_network_architectures.building_blocks.residual_encoders import (
    ResidualEncoder,
)
from dynamic_network_architectures.building_blocks.unet_decoder import UNetDecoder


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


class ResEncLUNet(nn.Module):
    """ResEncL encoder + UNetDecoder, 2-class (background / plaque) logits."""

    def __init__(
        self,
        num_classes: int = 2,
        num_input_channels: int = 1,
        n_conv_per_stage_decoder: int = 1,
        deep_supervision: bool = False,
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

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.decoder(self.encoder(x))
