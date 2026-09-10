"""Compatibility helpers for 5D (t, c, z, y, x) OME-Zarr stores."""

import math

from zarrnii import ZarrNii
from zarrnii.core import _derive_ngff_image


def load_near_isotropic(path, level, target_scale, **from_file_kwargs):
    """Load *path* on the grid nearest an isotropic ``target_scale``.

    ``ZarrNii.from_file(..., downsample_near_isotropic=True)`` cannot do this:
    it equalizes finer axes to the *coarsest* axis (whatever resolution that
    happens to be, not a chosen target) and refuses to act when more than one
    axis needs correcting. A pyramid that downsamples all three axes puts z at
    twice the target already at level 1, which no amount of further
    downsampling can undo -- the level itself has to change.

    Starting from the requested *level*, this walks to finer levels while any
    spatial axis is coarser than ``target_scale * sqrt(2)`` (past that point
    the next-finer level is strictly closer in log space), then block-average
    downsamples each remaining finer axis by the power of two nearest the
    target. The requested level is a cap on cost, not a promise: the result is
    never coarser than it, but may be finer on some axes.

    ``target_scale`` is in the units of the store's scale metadata (µm for
    SPIMprep outputs).
    """
    while True:
        znimg = ZarrNii.from_file(path, level=level, **from_file_kwargs)
        spatial = [d for d in ("z", "y", "x") if d in znimg.scale]
        scales = [float(znimg.scale[d]) for d in spatial]
        if level == 0 or all(s <= target_scale * math.sqrt(2) for s in scales):
            break
        level -= 1

    factors = [2 ** max(0, round(math.log2(target_scale / s))) for s in scales]
    if any(f > 1 for f in factors):
        znimg = znimg.downsample(factors=factors, spatial_dims=spatial)

    coarse = {
        d: s * f
        for d, s, f in zip(spatial, scales, factors)
        if s * f > target_scale * math.sqrt(2)
    }
    if coarse:
        print(
            f"WARNING: axes {coarse} are coarser than the {target_scale} target "
            f"even at level 0; proceeding on the finest grid available",
            flush=True,
        )
    return znimg


def drop_singleton_time(znimg):
    """Return *znimg* without its time axis, which must be singleton.

    Some acquisitions store a leading singleton t axis.
    ``apply_scaled_processing`` aligns the hires and lowres images by
    positional dimension index, so a 5D hires image against a 4D lowres
    bias field mis-slices the lowres array (map_coordinates raises
    "invalid shape for coordinate array"). Dropping the axis up front
    puts 5D stores on the same (c, z, y, x) path as every other dataset.
    """
    if "t" not in znimg.dims:
        return znimg
    t = znimg.dims.index("t")
    if znimg.data.shape[t] != 1:
        raise ValueError(
            f"non-singleton time axis (size {znimg.data.shape[t]}) is not supported"
        )
    index = tuple(0 if i == t else slice(None) for i in range(znimg.data.ndim))
    ngff = znimg.ngff_image
    new_ngff = _derive_ngff_image(
        ngff,
        data=znimg.data[index],
        dims=[d for d in ngff.dims if d != "t"],
        scale={k: v for k, v in ngff.scale.items() if k != "t"},
        translation={k: v for k, v in ngff.translation.items() if k != "t"},
    )
    return ZarrNii(
        ngff_image=new_ngff,
        axes_order=znimg.axes_order,
        xyz_orientation=znimg.xyz_orientation,
        _omero=znimg._omero,
    )
