"""Extract a user-defined bounding box crop from SPIM zarr data in template space.

Resamples the SPIM zarr at the requested input level to a user-specified bounding
box in template space at the desired isotropic resolution, saving as a NIfTI file.

The reference volume in template space is constructed from the bounding box centre
coordinates (RAS mm) and the target voxel size, then the floating SPIM zarr is
pull-resampled to that reference using the warp (subject →  template)
via zarrnii's block-wise interpolation (ZarrNii apply_transform)

A special case (use_brain_mask: true) uses the template brain mask to define the
bounding box extent, enabling whole-brain resampling at arbitrary resolution.

Bounding box config format
--------------------------
Each entry in the bbox YAML file (passed via --bbox_config) can have two forms:

Explicit coordinate form::

    - name: entorhinal
      x: 5.2          # centre in RAS mm (R-axis)
      y: -2.3         # centre in RAS mm (A-axis)
      z: -1.8         # centre in RAS mm (S-axis)
      units: mm       # optional; 'mm' (default) or 'um'
      size: [256, 256, 256]  # [nx, ny, nz] voxels along [R, A, S]
      resolution_um: 20      # isotropic output resolution in microns
      input_level: 0         # which zarr pyramid level to read from

Brain-mask form::

    - name: wholebrain
      use_brain_mask: true
      resolution_um: 100
      input_level: 2

Coordinate system notes
-----------------------
* ``x``, ``y``, ``z`` are in standard RAS convention: x=Right, y=Anterior, z=Superior.
* ``size`` = [nx, ny, nz] where nx/ny/nz are voxel counts along the R/A/S axes
  respectively.

"""

import dask.array as da
import nibabel as nib
import numpy as np
from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from zarrnii import DisplacementTransform, ZarrNii

bbox_config = snakemake.params.bbox_config
stain = snakemake.wildcards.stain
input_level = int(snakemake.wildcards.level)
zarrnii_kwargs = snakemake.params.zarrnii_kwargs

# Load floating SPIM zarr (ZYX axes_order = zarrnii default)
flo_znimg = ZarrNii.from_file(
    snakemake.input.spim,
    level=input_level,
    channel_labels=[stain],
    **zarrnii_kwargs,
)

# Build reference geometry (all coordinates in RAS mm order)
resolution_um = bbox_config["resolution_um"]
vox_size_mm = resolution_um / 1000.0

if bbox_config.get("use_brain_mask", False):
    # Derive bounding box from the template brain mask
    mask_nib = nib.load(snakemake.input.brain_mask)
    mask_affine = mask_nib.affine
    mask_data = mask_nib.get_fdata()

    # Find the voxel extent of non-zero mask voxels
    nz_indices = np.argwhere(mask_data > 0)
    min_vox = nz_indices.min(axis=0)
    max_vox = nz_indices.max(axis=0)

    # Convert both corners to RAS mm (handling possible negative affine strides)
    corners_mm = np.column_stack(
        [
            mask_affine[:3, :3] @ np.column_stack([min_vox, max_vox])
            + mask_affine[:3, 3:4]
        ]
    ).T  # shape (2, 3): each row is (R, A, S) mm

    origin_mm = corners_mm.min(axis=0)  # (R_start, A_start, S_start)
    max_mm = corners_mm.max(axis=0)
    extent_mm = max_mm - origin_mm

    nx = int(np.ceil(extent_mm[0] / vox_size_mm))  # voxels along R
    ny = int(np.ceil(extent_mm[1] / vox_size_mm))  # voxels along A
    nz = int(np.ceil(extent_mm[2] / vox_size_mm))  # voxels along S
else:
    # Explicit centre coordinates in template space
    units = bbox_config.get("units", "mm")
    cx = float(bbox_config["x"])
    cy = float(bbox_config["y"])
    cz = float(bbox_config["z"])
    if units == "um":
        cx, cy, cz = cx / 1000.0, cy / 1000.0, cz / 1000.0
    elif units != "mm":
        raise ValueError(
            f"Unsupported units '{units}' in bbox_config; use 'mm' or 'um'"
        )

    size = bbox_config.get("size", [256, 256, 256])
    nx, ny, nz = int(size[0]), int(size[1]), int(size[2])

    origin_mm = np.array(
        [
            cx - nx * vox_size_mm / 2.0,  # R_start
            cy - ny * vox_size_mm / 2.0,  # A_start
            cz - nz * vox_size_mm / 2.0,  # S_start
        ]
    )

ref_darr = da.zeros((1, nz, ny, nx), dtype=np.float32, chunks=(1, 64, 64, 64))
ref_znimg = ZarrNii.from_darr(
    ref_darr,
    axes_order="ZYX",
    orientation="RAS",
    spacing=(vox_size_mm, vox_size_mm, vox_size_mm),
    origin=(float(origin_mm[2]), float(origin_mm[1]), float(origin_mm[0])),
    axes_units={"x": "millimeter", "y": "millimeter", "z": "millimeter"},
    channel_labels=[stain],
)

# Load composite warp (template → subject)
composite_transform = DisplacementTransform.from_nifti(snakemake.input.xfm_composite)


out_znimg = flo_znimg.apply_transform(composite_transform, ref_znimg=ref_znimg)

with ProgressBar():
    out_znimg.to_nifti(snakemake.output.nii)
