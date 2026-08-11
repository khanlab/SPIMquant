"""Threshold the LANTERN vote map and upsample it to the segmentation level.

The ensemble runs at ``plaque_level`` because it is scale-sensitive and was built
for that grid; the segmentation level output is therefore the upsampled
prediction, NOT native inference at that level.

Upsampling is exact replication (``da.repeat``) rather than interpolation, and
the per-axis factors are derived from the shape ratio so this stays correct when
the pyramid downsamples only x and y.

Output is valued 0/100 rather than 0/1, matching ``gmmthresh`` / ``multiotsu`` /
``threshold`` in segmentation.smk, so that mean-pool downsampling in ``fieldfrac``
yields a percentage directly.
"""

import dask.array as da
import numpy as np
from dask.diagnostics import ProgressBar
from zarrnii import ZarrNii


def upsample_nearest(arr, target_shape):
    """Replicate voxels to `target_shape`, then trim/pad to land on it exactly."""
    if len(target_shape) != arr.ndim:
        raise ValueError(f"rank mismatch: {arr.shape} vs target {target_shape}")

    out = arr
    for axis, (have, want) in enumerate(zip(arr.shape, target_shape)):
        factor = max(1, int(round(want / have)))
        if factor > 1:
            out = da.repeat(out, factor, axis=axis)

    out = out[tuple(slice(0, min(h, w)) for h, w in zip(out.shape, target_shape))]
    pad = [(0, w - h) for h, w in zip(out.shape, target_shape)]
    if any(after for _, after in pad):
        out = da.pad(out, pad, mode="edge")
    return out


def main():
    n_folds = int(snakemake.params.n_folds)
    vote_threshold = int(snakemake.params.vote_threshold)

    if not 1 <= vote_threshold <= n_folds:
        raise ValueError(
            f"plaque_vote_threshold={vote_threshold} is out of range for an "
            f"{n_folds}-fold ensemble; it must be between 1 and {n_folds}"
        )

    probseg = ZarrNii.from_file(snakemake.input.probseg, level=0)
    ref = ZarrNii.from_file(
        snakemake.input.spim,
        level=int(snakemake.params.target_level),
        channel_labels=[snakemake.wildcards.stain],
        **snakemake.params.zarrnii_kwargs,
    )

    # probseg holds votes/n_folds. Compare against the MIDPOINT between adjacent
    # attainable values rather than the exact fraction: 2.5/5 = 0.5 sits safely
    # between 0.4 and 0.6, so no float representation error can flip a voxel at the
    # boundary the way `>= 3/5` could.
    binary = probseg.data >= (vote_threshold - 0.5) / n_folds

    upsampled = upsample_nearest(binary, ref.data.shape)

    znimg_mask = ref.copy()
    znimg_mask.data = (upsampled * 100).astype(np.uint8)

    print(
        f"votes>={vote_threshold} of {n_folds} | "
        f"{probseg.data.shape} -> {ref.data.shape}",
        flush=True,
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
