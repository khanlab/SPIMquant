"""Fibre orientation maps from moments (rule fibre_metrics).

Reads the per-super-voxel moments written by fibre_moments.py and writes:

* one scalar NIfTI per metric (fibre density, mean Frangi response, mean
  scale, dispersion indices, orientation-tensor FA, ODF GFA); orientation
  metrics are NaN in blocks with too little fibre (``fibre_min_voxels``),
* the principal fibre direction (3 volumes, world x/y/z) and its
  direction-encoded colour map (|direction| x FA),
* the ODF as even-order real SH coefficients in the MRtrix3 convention
  (view with ``mrview -odf.load_sh``).
"""

import numpy as np

SCALAR_METRICS = (
    "fibredensity",
    "frangi",
    "scale",
    "odi",
    "odipri",
    "odisec",
    "odianis",
    "fa",
    "gfa",
)


# float64 working memory per slab of moments read at a time (bytes)
_SLAB_BYTES = 5e8


def compute_maps(moments, odf=True, min_fibre=0.0, missing=0.0):
    """All maps from a moment array (n_moments, z, y, x) as numpy.

    Returns a dict of ``SCALAR_METRICS`` (3D), ``dir`` and ``dec``
    (3, z, y, x) and, if ``odf``, ``odf`` (n_coeffs, z, y, x).  Orientation
    metrics are ``missing`` where a block has at most ``min_fibre`` weighted
    fibre voxels (see ``moments_to_metrics``); the DEC map is 0 there.  The
    ODF is not masked.
    """
    from fibre_orientation import moments_to_metrics, odf_coefficients

    m = moments_to_metrics(moments, min_fibre=min_fibre, missing=missing)
    out = {k: m[k] for k in SCALAR_METRICS}
    direction = np.stack([m["dirx"], m["diry"], m["dirz"]])
    out["dir"] = direction
    out["dec"] = np.nan_to_num(np.abs(direction) * m["fa"][None]).astype(np.float32)
    if odf:
        out["odf"] = odf_coefficients(moments).astype(np.float32)
    return out


def z_slabs(darr, slab_bytes=_SLAB_BYTES):
    """Yield ``(z0, z1, moments)`` float64 slabs of a (c, z, y, x) array."""
    c, nz, ny, nx = darr.shape
    step = max(1, int(slab_bytes // (8 * c * ny * nx)))
    for z0 in range(0, nz, step):
        z1 = min(z0 + step, nz)
        yield z0, z1, np.asarray(darr[:, z0:z1], dtype=np.float64)


def write_nifti(template, data, path):
    """Write ``data`` (c, z, y, x) or (z, y, x) on the grid of ``template``."""
    import dask.array as da

    zn = template.copy()
    arr = np.asarray(data, dtype=np.float32)
    if arr.ndim == 3:
        arr = arr[None]
    zn.data = da.from_array(arr)
    zn.to_nifti(path)


def main():
    import dask
    from fibre_orientation import SH_START, odf_coefficients
    from zarrnii import ZarrNii

    out = snakemake.output
    with dask.config.set(scheduler="threads", num_workers=snakemake.threads):
        znimg = ZarrNii.from_file(snakemake.input.moments)
        darr = znimg.data
        _, nz, ny, nx = darr.shape

        # scalar maps, direction and DEC (float32, filled slab by slab)
        names = list(SCALAR_METRICS)
        maps = np.zeros((len(names) + 6, nz, ny, nx), dtype=np.float32)
        for z0, z1, m in z_slabs(darr):
            res = compute_maps(
                m,
                odf=False,
                min_fibre=float(snakemake.params.min_fibre),
                missing=np.nan,
            )
            for i, name in enumerate(names):
                maps[i, z0:z1] = res[name]
            maps[len(names) : len(names) + 3, z0:z1] = res["dir"]
            maps[len(names) + 3 :, z0:z1] = res["dec"]
        for i, name in enumerate(names):
            write_nifti(znimg, maps[i], getattr(out, name))
        write_nifti(znimg, maps[len(names) : len(names) + 3], out.dir)
        write_nifti(znimg, maps[len(names) + 3 :], out.dec)
        del maps

        # ODF (SH coefficients)
        odf = np.zeros((darr.shape[0] - SH_START, nz, ny, nx), dtype=np.float32)
        for z0, z1, m in z_slabs(darr):
            odf[:, z0:z1] = odf_coefficients(m)
        write_nifti(znimg, odf, out.odf)
        del odf


if __name__ == "__main__":
    main()
