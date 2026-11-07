"""Brain-masked intensities without bias-field correction.

``--correction_method none``: the raw image, zeroed outside the brain mask.
``--correction_method log``:  ``log1p`` of the raw image, zeroed outside the
brain mask.

Bias-field correction (n4, gaussian) fits a smooth field at a coarse level,
upsamples it to full resolution and divides it out, which is the heaviest
step before segmentation.  Methods that respond to local contrast, such as
DoG blobs, do not need it; in log space a multiplicative bias becomes an
additive offset that the DoG background term removes.

The brain mask is still applied, because downstream steps treat zero as
"outside the brain" (e.g. the mask-normalised blurs in dogthresh and the
instance splitting).  It is upsampled to each block by linear interpolation
and thresholded at 0.5, as for n4, so the image is read and written once with
no full-resolution intermediates.
"""

import numpy as np


def _axis_weights(start, size, n_full, n_low):
    """Linear-interpolation indices and weights along one axis, for full-res
    voxels start..start+size of n_full, on a grid of n_low voxels covering
    the same extent (voxel centres aligned)."""
    pos = (np.arange(start, start + size) + 0.5) * (n_low / n_full) - 0.5
    pos = np.clip(pos, 0, n_low - 1)
    i0 = np.floor(pos).astype(int)
    i1 = np.minimum(i0 + 1, n_low - 1)
    return i0, i1, (pos - i0).astype(np.float32)


def upsample_mask_block(lowres_mask, block_shape, block_start, full_shape):
    """Linearly interpolated low-resolution mask for one full-resolution block,
    thresholded at 0.5 (as n4_biasfield's mask is).  Separable, and reads only
    the low-resolution voxels the block needs."""
    vals = np.asarray(lowres_mask, dtype=np.float32)
    for axis, (start, size, n_full) in enumerate(
        zip(block_start, block_shape, full_shape)
    ):
        i0, i1, w = _axis_weights(start, size, n_full, lowres_mask.shape[axis])
        shape = [1] * vals.ndim
        shape[axis] = -1
        w = w.reshape(shape)
        vals = (1 - w) * np.take(vals, i0, axis=axis) + w * np.take(vals, i1, axis=axis)
    return vals > 0.5


def mask_block(block, lowres_mask, full_shape, method, block_info=None):
    """One block of the output: masked raw (none) or masked log1p (log)."""
    start = [loc[0] for loc in block_info[0]["array-location"]]
    inside = upsample_mask_block(lowres_mask, block.shape, start, full_shape)
    if method == "log":
        vals = np.log1p(np.maximum(block, 0).astype(np.float32))
        return np.where(inside, vals, np.float32(0))
    return np.where(inside, block, np.zeros((), dtype=block.dtype))


def check_mask_shape(mask_shape, full_shape, factor, tolerance=2):
    """The mask must be the image downsampled by ``factor`` on every axis."""
    expected = [n / factor for n in full_shape]
    if any(abs(m - e) > tolerance for m, e in zip(mask_shape, expected)):
        raise ValueError(
            f"Brain mask shape {tuple(mask_shape)} does not match the image "
            f"{tuple(full_shape)} downsampled by {factor:g} "
            f"(expected about {tuple(round(e) for e in expected)})"
        )


if __name__ == "__main__":
    import dask.array as da
    from dask_setup import get_dask_client
    from zarrnii import ZarrNii
    from zarrnii_compat import drop_singleton_time

    method = snakemake.wildcards.method
    if method not in ("none", "log"):
        raise ValueError(f"Unknown correction method {method!r}")
    is_imaris = str(snakemake.input.spim).lower().endswith(".ims")

    with get_dask_client(
        snakemake.config["dask_scheduler"],
        snakemake.threads,
        threads_per_worker=16 if is_imaris else 2,
    ):
        hires_level = int(snakemake.wildcards.level)
        mask_level = int(snakemake.params.mask_level)

        znimg = ZarrNii.from_file(
            snakemake.input.spim,
            channel_labels=[snakemake.wildcards.stain],
            level=hires_level,
            downsample_near_isotropic=True,
            chunks=(256, 256, 256) if is_imaris else None,
            **snakemake.params.zarrnii_kwargs,
        )
        znimg = drop_singleton_time(znimg)
        if znimg.data.ndim != 4 or znimg.data.shape[0] != 1:
            raise ValueError(
                f"Expected a single-channel (c, z, y, x) image, got {znimg.data.shape}"
            )
        img = znimg.data[0]

        lowres_mask = (
            np.asarray(
                ZarrNii.from_nifti(snakemake.input.mask, axes_order="ZYX").data
            ).squeeze()
            > 0
        )
        check_mask_shape(lowres_mask.shape, img.shape, 2 ** (mask_level - hires_level))
        print(
            f"{method}: image {img.shape} {img.dtype}, brain mask {lowres_mask.shape}, "
            f"{100 * lowres_mask.mean():.1f}% inside"
        )

        out_dtype = np.float32 if method == "log" else img.dtype
        corrected = da.map_blocks(
            mask_block,
            img,
            lowres_mask=lowres_mask,
            full_shape=img.shape,
            method=method,
            dtype=out_dtype,
            meta=np.array((), dtype=out_dtype),
        )

        znimg_out = znimg.copy()
        znimg_out.data = corrected[np.newaxis]
        znimg_out.to_ome_zarr(
            snakemake.output.corrected,
            match_scale_factors_from=snakemake.input.spim,
            **snakemake.config["zarrnii_out_kwargs"],
        )
