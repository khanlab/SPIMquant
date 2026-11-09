"""Fibre orientation maps from moments (rule fibre_metrics).

Reads the per-super-voxel moments written by fibre_moments.py and writes:

* one scalar NIfTI per metric (fibre density, mean Frangi response, mean
  scale, dispersion indices, orientation-tensor FA, ODF GFA),
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


def compute_maps(moments):
    """All maps from a moment array (n_moments, z, y, x) as numpy.

    Returns a dict of ``SCALAR_METRICS`` (3D), ``dir`` and ``dec``
    (3, z, y, x) and ``odf`` (n_coeffs, z, y, x).
    """
    from fibre_orientation import moments_to_metrics, odf_coefficients

    m = moments_to_metrics(moments)
    out = {k: m[k] for k in SCALAR_METRICS}
    direction = np.stack([m["dirx"], m["diry"], m["dirz"]])
    out["dir"] = direction
    out["dec"] = (np.abs(direction) * m["fa"][None]).astype(np.float32)
    out["odf"] = odf_coefficients(moments).astype(np.float32)
    return out


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
    from zarrnii import ZarrNii

    with dask.config.set(scheduler="threads", num_workers=snakemake.threads):
        znimg = ZarrNii.from_file(snakemake.input.moments)
        moments = np.asarray(znimg.data.compute(), dtype=np.float64)

    maps = compute_maps(moments)
    single = znimg.select_channels([0]) if hasattr(znimg, "select_channels") else None
    for name in SCALAR_METRICS:
        write_nifti(single or znimg, maps[name], getattr(snakemake.output, name))
    write_nifti(znimg, maps["dir"], snakemake.output.dir)
    write_nifti(znimg, maps["dec"], snakemake.output.dec)
    write_nifti(znimg, maps["odf"], snakemake.output.odf)


if __name__ == "__main__":
    main()
