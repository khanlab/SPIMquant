"""Resample SPIM OME-Zarr data into template space at a given resolution.

The output grid has isotropic voxels of exactly the requested resolution
(``{res}`` in um, with p for a decimal point, e.g. 2p5 for 2.5 um) and covers
the registration template's field of view, extending past it by less than a
voxel where the resolution does not divide the template's extent. Each output
voxel is pulled from the SPIM image (at the requested pyramid level) through
the subject-to-template composite warp, using zarrnii's block-wise resampler,
and the result is written as an OME-Zarr pyramid.

The input is either the raw SPIM channel (``desc-raw``) or its intensity
corrected image at the same level (``desc-corrected<method>``).
"""

import math
import shutil
import tempfile
import time
from pathlib import Path

import dask.array as da
import nibabel as nib
import numpy as np
from dask_setup import get_dask_client
from zarrnii import DisplacementTransform, ZarrNii


def template_grid(path, res_mm, chunks=128):
    """Lazy ZYX grid with res_mm voxels over a NIfTI template's field of view.

    The grid has the template's orientation and starts at the outer corner of
    its first voxel; each axis has as many res_mm voxels as cover the
    template's extent, so the grid can extend past the template by less than
    a voxel. Built from the NIfTI affine rather than
    ZarrNii.from_nifti(axes_order="ZYX"), which in zarrnii <= 0.22 swaps the x
    and z offsets (and voxel sizes).
    """
    img = nib.load(path)
    aff = img.affine
    zooms = img.header.get_zooms()[:3]
    shape, origin = [], []
    for j in range(3):  # NIfTI voxel axes i, j, k
        row = int(np.argmax(np.abs(aff[:3, j])))  # world axis of voxel axis j
        sign = np.sign(aff[row, j])
        extent = img.shape[j] * float(zooms[j])
        # an exact multiple (up to float32 voxel sizes in NIfTI headers) does
        # not get an extra voxel
        n = extent / res_mm
        shape.append(max(1, round(n) if abs(n - round(n)) < 1e-3 else math.ceil(n)))
        # centre of the first voxel: template corner + half a new voxel
        origin.append(float(aff[row, 3]) + sign * (res_mm - float(zooms[j])) / 2)
    nx, ny, nz = shape
    return ZarrNii.from_darr(
        da.zeros((1, nz, ny, nx), dtype=np.float32, chunks=(1, chunks, chunks, chunks)),
        axes_order="ZYX",
        orientation="".join(nib.aff2axcodes(aff)),
        spacing=(res_mm, res_mm, res_mm),
        origin=tuple(origin[::-1]),
        axes_units={"x": "millimeter", "y": "millimeter", "z": "millimeter"},
    )


def coarse_chunks(znimg, target_mm=0.4):
    """Chunks for *znimg* made of whole native chunks, about target_mm per side.

    Each output block reads the input chunks it overlaps, so tiny native
    chunks (e.g. 32x128x128 in Imaris files) mean hundreds of thousands of
    input tasks and many per output block.
    """
    spatial = znimg.data.chunksize[-3:]
    factors = [
        max(1, round(target_mm / (c * float(sp))))
        for c, sp in zip(spatial, znimg.get_zooms())
    ]
    lead = (1,) * (znimg.data.ndim - 3)
    return lead + tuple(c * f for c, f in zip(spatial, factors))


def pyramid_levels(shape, coarsest=128):
    """Number of pyramid levels (incl. level 0), halving until every axis is at
    most ``coarsest`` voxels, as Imaris does (e.g. 7 for 3200x5280x4560)."""
    return 1 + max(0, math.ceil(math.log2(max(shape) / coarsest)))


def main():
    wildcards = snakemake.wildcards
    params = snakemake.params
    res_mm = float(wildcards.res.replace("p", ".")) / 1000

    with get_dask_client(
        snakemake.config["dask_scheduler"],
        snakemake.threads,
        mem_mb=snakemake.resources.mem_mb,
        reserve_mb=4000,
    ):
        if params.corrected:
            # already single channel, with orientation applied when it was written
            flo = ZarrNii.from_file(snakemake.input.spim)
        else:
            load = dict(
                level=int(wildcards.level),
                channel_labels=[wildcards.stain],
                **params.zarrnii_kwargs,
            )
            flo = ZarrNii.from_file(snakemake.input.spim, **load)
            chunks = coarse_chunks(flo)
            if chunks != flo.data.chunksize:
                print(f"reading input chunks {flo.data.chunksize} as {chunks}")
                flo = ZarrNii.from_file(snakemake.input.spim, chunks=chunks, **load)

        grid = template_grid(snakemake.input.ref, res_mm, chunks=params.chunks)
        # subject -> template warp, defined on the template grid: maps template
        # points to subject points, as used by antsApplyTransforms -r template.
        # Saved to zarr so that each output block reads only the part it needs,
        # instead of the whole field going into the dask graph.
        warp_dir = Path(
            tempfile.mkdtemp(suffix="_warp.zarr", dir=snakemake.resources.tmpdir)
        )
        warp = DisplacementTransform.from_nifti(snakemake.input.xfm_composite).to_zarr(
            warp_dir
        )

        t0 = time.time()
        out = flo.apply_transform(
            warp,
            ref_znimg=grid,
            chunks=params.chunks,
            io=params.io,
        )
        print(
            f"input {flo.shape} {flo.data.dtype} (spacing {flo.get_zooms()} mm) -> "
            f"output {out.shape} {out.data.dtype} (spacing {out.get_zooms()} mm), "
            f"{out.data.numblocks} blocks, planned in {time.time() - t0:.1f}s",
            flush=True,
        )

        t0 = time.time()
        try:
            out.to_ome_zarr(
                snakemake.output.spim,
                max_layer=pyramid_levels(out.shape[-3:]),
                **snakemake.config["zarrnii_out_kwargs"],
            )
        finally:
            shutil.rmtree(warp_dir, ignore_errors=True)
        print(f"resampled and written in {time.time() - t0:.1f}s", flush=True)


# the distributed scheduler's worker processes import this file: only the
# snakemake job itself should run it
if __name__ == "__main__":
    main()
