"""Resample SPIM OME-Zarr data into template space, on a zoomed template grid.

The output grid covers the registration template's field of view with voxels
``zoom`` times smaller than the template's (e.g. zoom 10 on a 25 um template
gives 2.5 um voxels). Each output voxel is pulled from the SPIM image (at the
requested pyramid level) through the subject-to-template composite warp, using
zarrnii's block-wise resampler, and the result is written as an OME-Zarr
pyramid.

The input is either the raw SPIM channel (``desc-raw``) or its intensity
corrected image at the same level (``desc-corrected<method>``).
"""

import shutil
import tempfile
import time
from pathlib import Path

import dask.array as da
import nibabel as nib
import numpy as np
from dask_setup import get_dask_client
from zarrnii import DisplacementTransform, ZarrNii


def template_grid(path, chunks=64):
    """Lazy ZYX reference grid with the geometry of a NIfTI image.

    Built from the NIfTI affine rather than ZarrNii.from_nifti(axes_order="ZYX"),
    which in zarrnii <= 0.22 swaps the x and z offsets (and voxel sizes).
    """
    img = nib.load(path)
    aff = img.affine
    zooms = img.header.get_zooms()[:3]
    # for each voxel axis (i, j, k): its world axis and the offset along it
    rows = [int(np.argmax(np.abs(aff[:3, j]))) for j in range(3)]
    offsets = [float(aff[rows[j], 3]) for j in range(3)]
    nx, ny, nz = img.shape[:3]
    return ZarrNii.from_darr(
        da.zeros((1, nz, ny, nx), dtype=np.float32, chunks=(1, chunks, chunks, chunks)),
        axes_order="ZYX",
        orientation="".join(nib.aff2axcodes(aff)),
        spacing=tuple(float(z) for z in zooms[::-1]),
        origin=tuple(offsets[::-1]),
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


def main():
    wildcards = snakemake.wildcards
    params = snakemake.params
    zoom = float(wildcards.zoom.replace("p", "."))

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

        template = template_grid(snakemake.input.ref)
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
            ref_znimg=template,
            ref_zoom=zoom,
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
                snakemake.output.spim, **snakemake.config["zarrnii_out_kwargs"]
            )
        finally:
            shutil.rmtree(warp_dir, ignore_errors=True)
        print(f"resampled and written in {time.time() - t0:.1f}s", flush=True)


# the distributed scheduler's worker processes import this file: only the
# snakemake job itself should run it
if __name__ == "__main__":
    main()
