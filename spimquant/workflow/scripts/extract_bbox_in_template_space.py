"""Extract a user-defined bounding box crop from SPIM zarr data in template space.

Resamples the SPIM zarr at the requested input level to a user-specified bounding
box in template space at the desired isotropic resolution, saving as a NIfTI file.

The reference volume in template space is constructed from the bounding box centre
coordinates (RAS mm) and the target voxel size, then the floating SPIM zarr is
pull-resampled to that reference using the composite inverse warp (template → subject)
via zarrnii's block-wise interpolation.

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

Transform pipeline
------------------
The block-wise resampling pipeline (ZYX voxel ordering throughout):

1. ``ref ZYX vox → RAS mm`` via the reference ZYX affine from zarrnii
   (``get_affine_matrix("ZYX")`` outputs (R, A, S) mm for both ZYX and XYZ order)
2. ``composite_inv`` displacement: RAS template mm → RAS subject mm
3. ``flo ZYX aff⁻¹``: RAS subject mm → flo ZYX vox

No coordinate permutation is needed because zarrnii's ``get_affine_matrix("ZYX")``
already maps ZYX voxel indices to (R, A, S) mm, which is what
``DisplacementTransform`` expects.
"""

import dask.array as da
import nibabel as nib
import numpy as np
from dask.diagnostics import ProgressBar
from dask_setup import get_dask_client
from zarrnii import AffineTransform, DisplacementTransform, ZarrNii
from zarrnii.core import interp_by_block

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

# Load composite inverse warp (template → subject)
disp_transform = DisplacementTransform.from_nifti(snakemake.input.xfm_composite_inv)

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

# Build reference ZarrNii using from_darr (ZYX axes_order, RAS orientation).
#
# from_darr(axes_order="ZYX", origin=(oz, oy, ox)) sets:
#   translation = {"z": oz, "y": oy, "x": ox}
# get_affine_matrix("ZYX") with RAS orientation uses reversed axcodes "SAR":
#   R_mm = X_idx * sx + ox   (x ↔ R)
#   A_mm = Y_idx * sy + oy   (y ↔ A)
#   S_mm = Z_idx * sz + oz   (z ↔ S)
#
# So origin should be (S_start, A_start, R_start) = (origin_mm[2], origin_mm[1], origin_mm[0])
# and the reference array shape is (C=1, Z=nz, Y=ny, X=nx).

ref_darr = da.zeros((1, nz, ny, nx), dtype=np.float32, chunks=(1, 64, 64, 64))
ref_znimg = ZarrNii.from_darr(
    ref_darr,
    axes_order="ZYX",
    orientation="RAS",
    spacing=(vox_size_mm, vox_size_mm, vox_size_mm),
    origin=(float(origin_mm[2]), float(origin_mm[1]), float(origin_mm[0])),
    channel_labels=[stain],
)

# Build transform chain (ZYX voxel space throughout):
#
# zarrnii's get_affine_matrix("ZYX") maps (Z_idx, Y_idx, X_idx) → (R_mm, A_mm, S_mm),
# i.e. the output is in standard RAS physical order — exactly what DisplacementTransform
# expects (it operates in RAS mm and uses the NIfTI affine to convert to disp-field
# voxels).
#
#   Step 1: ref ZYX vox  ──(ref_zyx_aff)──►  RAS mm  (template space)
#   Step 2: RAS mm       ──(disp_transform)──►  RAS mm  (subject space)
#   Step 3: RAS mm       ──(inv flo_zyx_aff)──►  flo ZYX vox

ref_zyx_aff = ref_znimg.get_affine_matrix("ZYX")
flo_zyx_aff = flo_znimg.get_affine_matrix("ZYX")

transforms = [
    AffineTransform.from_array(ref_zyx_aff),  # ref ZYX vox → RAS mm
    disp_transform,  # RAS template → RAS subject
    AffineTransform.from_array(np.linalg.inv(flo_zyx_aff)),  # RAS mm → flo ZYX vox
]

# Prefer direct zarr access to avoid nested dask compute() calls
store_info = flo_znimg.get_zarr_store_info()
if store_info is not None:
    flo_store_path = str(store_info["store_path"])
    flo_array_shape = store_info["array_shape"]
    flo_dataset_path = store_info["dataset_path"]
    flo_znimg_arg = None
else:
    # Fall back to in-memory path (legacy behaviour)
    flo_store_path = None
    flo_array_shape = None
    flo_dataset_path = "0"
    flo_znimg_arg = flo_znimg

with get_dask_client("threads", snakemake.threads):
    resampled = da.map_blocks(
        interp_by_block,
        ref_znimg.darr,
        dtype=np.float32,
        transforms=transforms,
        flo_store_path=flo_store_path,
        flo_array_shape=flo_array_shape,
        flo_dataset_path=flo_dataset_path,
        flo_znimg=flo_znimg_arg,
    )

    # Wrap result in ZarrNii with the same geometry as the reference
    result_znimg = ZarrNii.from_darr(
        resampled,
        axes_order=ref_znimg.axes_order,
        orientation=ref_znimg.xyz_orientation,
        spacing=(vox_size_mm, vox_size_mm, vox_size_mm),
        origin=(float(origin_mm[2]), float(origin_mm[1]), float(origin_mm[0])),
        channel_labels=[stain],
    )

    with ProgressBar():
        result_znimg.to_nifti(snakemake.output.nii)
