"""Compatibility helpers for 5D (t, c, z, y, x) OME-Zarr stores and Imaris files."""

import math

from zarrnii import ZarrNii
from zarrnii.core import _derive_ngff_image

# Millimetres per unit, for every OME-Zarr space unit zarrnii may report.
_MM_PER_UNIT = {
    "nanometer": 1e-6,
    "nm": 1e-6,
    "micrometer": 1e-3,
    "micron": 1e-3,
    "um": 1e-3,
    "millimeter": 1.0,
    "mm": 1.0,
    "centimeter": 10.0,
    "cm": 10.0,
    "meter": 1e3,
    "m": 1e3,
}


def _um_to_axis_units(znimg, value_um, axis):
    """Convert *value_um* (micrometres) into the units *znimg* stores *axis* in.

    zarrnii normalises OME-Zarr and Imaris scales to millimetres and records the
    unit per axis in ``ngff_image.axes_units``; an axis with no recorded unit is
    taken to be in millimetres for the same reason. Interpreting a micrometre
    target as if it were in the store's units is off by 1000x, and downsamples a
    whole brain to a handful of voxels.
    """
    units = getattr(znimg.ngff_image, "axes_units", None) or {}
    unit = (units.get(axis) or "millimeter").lower()
    if unit not in _MM_PER_UNIT:
        raise ValueError(f"unrecognised unit {unit!r} on axis {axis!r}")
    return value_um * 1e-3 / _MM_PER_UNIT[unit]


def axis_scale_um(znimg, axis):
    """The voxel size of *axis* in micrometres, whatever units the store uses."""
    return float(znimg.scale[axis]) / _um_to_axis_units(znimg, 1.0, axis)


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
    downsamples each axis that the walk made finer than it needs to be, by
    the power of two nearest the target. The requested level is a cap on cost
    and on coarseness: no axis is ever made coarser than it is at that level,
    so a level that is already near-isotropic at the target is returned as
    is, and only axes that the walk over-refined are downsampled back.

    ``target_scale`` is in micrometres regardless of the store's own units
    (zarrnii reports millimetres for both OME-Zarr and Imaris); it is converted
    per axis with ``_um_to_axis_units``.
    """
    cap = None
    while True:
        znimg = ZarrNii.from_file(path, level=level, **from_file_kwargs)
        spatial = [d for d in ("z", "y", "x") if d in znimg.scale]
        scales = [float(znimg.scale[d]) for d in spatial]
        targets = [_um_to_axis_units(znimg, target_scale, d) for d in spatial]
        if cap is None:
            cap = scales  # the requested level's own voxel size, per axis
        if level == 0 or all(
            s <= t * math.sqrt(2) for s, t in zip(scales, targets)
        ):
            break
        level -= 1

    factors = []
    for s, t, c_ in zip(scales, targets, cap):
        f = 2 ** max(0, round(math.log2(t / s)))
        while f > 1 and s * f > c_ * 1.001:
            f //= 2
        factors.append(f)
    if any(f > 1 for f in factors):
        znimg = znimg.downsample(factors=factors, spatial_dims=spatial)

    coarse = {
        d: s * f
        for d, s, f, t in zip(spatial, scales, factors, targets)
        if s * f > t * math.sqrt(2)
    }
    if coarse:
        print(
            f"WARNING: axes {coarse} (store units) are coarser than the "
            f"{target_scale} um target even at level 0; proceeding on the "
            f"finest grid available",
            flush=True,
        )
    return znimg


def open_inference_array(path, level, target_scale, **from_file_kwargs):
    """The lazy (c, z, y, x) array ``load_near_isotropic`` would run on.

    A singleton leading time axis, if the store has one, is dropped.  This is
    what reader processes call to open the input themselves, so that only the
    arguments (not a dask graph) have to be sent to them.
    """
    znimg = load_near_isotropic(path, level, target_scale, **from_file_kwargs)
    data = znimg.data
    if data.ndim == 5 and data.shape[0] == 1:
        data = data[0]
    return data


def with_channels(znimg, data, channel_labels):
    """A ZarrNii like *znimg* but carrying *data*, whose channel axis holds
    ``channel_labels`` (any number of channels, e.g. one probability map per
    predicted class).  Dims, scale and translation are kept; only the omero
    channel metadata is rebuilt so that the store advertises the new labels.
    """
    from zarrnii.core import make_omero

    ngff = znimg.ngff_image
    if data.ndim != len(ngff.dims):
        raise ValueError(
            f"data has {data.ndim} dims, image has {len(ngff.dims)} ({ngff.dims})"
        )
    new_ngff = _derive_ngff_image(ngff, data=data)
    return ZarrNii(
        ngff_image=new_ngff,
        axes_order=znimg.axes_order,
        xyz_orientation=znimg.xyz_orientation,
        _omero=make_omero(channel_labels=list(channel_labels)),
    )


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
