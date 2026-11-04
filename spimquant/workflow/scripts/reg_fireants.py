"""Register a moving image to a fixed image with FireANTs (GPU).

Writes ANTs-compatible forward and inverse displacement fields, which already
include the affine, so the affine file written for the shared compose step is
the identity. The warped images are the intensity-normalized images, for QC.
"""

import numpy as np
import SimpleITK as sitk

try:
    import torch
    from fireants.io.image import BatchedImages, Image
    from fireants.registration.affine import AffineRegistration
    from fireants.registration.greedy import GreedyRegistration
except ImportError as err:
    raise SystemExit(
        "FireANTs is not installed: run SPIMquant from the gpu environment, "
        "e.g. pixi run -e gpu spimquant ..."
    ) from err

p = snakemake.params
preset = p.preset
device = "cuda" if torch.cuda.is_available() else "cpu"


def load(path):
    """Load an image, clipped to robust percentiles and rescaled to [0, 1]."""
    image = Image.load_file(path, device=device)
    data = image.array.data
    values = data[data > 0]
    low, high = np.percentile(values.cpu().numpy(), preset["clip_percentiles"])
    image.array = ((data.clamp(low, high) - low) / max(high - low, 1e-6)).to(data.dtype)
    return image


def save_like(array, reference, path):
    out = sitk.GetImageFromArray(array[0, 0].detach().cpu().numpy())
    out.CopyInformation(reference.itk_image)
    sitk.WriteImage(out, path)


fixed_image = load(snakemake.input.fixed)
moving_image = load(snakemake.input.moving)
fixed = BatchedImages([fixed_image])
moving = BatchedImages([moving_image])


def center_alignment(fixed_image, moving_image):
    """Affine translating the fixed image's center to the moving image's center."""

    def center(image):
        itk = image.itk_image
        index = [(size - 1) / 2 for size in itk.GetSize()]
        return np.array(itk.TransformContinuousIndexToPhysicalPoint(index))

    matrix = torch.eye(4, device=device, dtype=torch.float32)[None]
    matrix[0, :3, 3] = torch.as_tensor(
        center(moving_image) - center(fixed_image), device=device, dtype=torch.float32
    )
    return matrix


init = center_alignment(fixed_image, moving_image) if p.init_centers else None

if preset["affine"]:
    affine = AffineRegistration(
        scales=preset["scales"],
        iterations=preset["affine_iters"],
        fixed_images=fixed,
        moving_images=moving,
        loss_type=preset["affine_loss"],
        init_rigid=init,
        progress_bar=False,
    )
    affine.optimize()
    init = affine.get_affine_matrix()
    with torch.no_grad():
        moved_affine = affine.evaluate(fixed, moving)
else:
    moved_affine = None

deform = GreedyRegistration(
    scales=preset["scales"],
    iterations=preset["deform_iters"],
    fixed_images=fixed,
    moving_images=moving,
    init_affine=init,
    loss_type=preset["deform_loss"],
    cc_kernel_size=preset["cc_kernel_size"],
    optimizer_lr=preset["optimizer_lr"],
    smooth_warp_sigma=preset["smooth_warp_sigma"],
    smooth_grad_sigma=preset["smooth_grad_sigma"],
    progress_bar=False,
)
deform.optimize()

with torch.no_grad():
    save_like(deform.evaluate(fixed, moving), fixed_image, snakemake.output.warped)
    if moved_affine is None:
        moved_affine = torch.nn.functional.grid_sample(
            moving(),
            torch.nn.functional.affine_grid(
                (moving.get_phy2torch() @ fixed.get_torch2phy())[:, :3],
                fixed.shape,
                align_corners=True,
            ),
            align_corners=True,
        )
    save_like(moved_affine, fixed_image, snakemake.output.warped_affine)

deform.save_as_ants_transforms([snakemake.output.warp])
deform.save_as_ants_transforms([snakemake.output.invwarp], save_inverse=True)

with open(snakemake.output.xfm_itk, "w") as f:
    f.write(
        "#Insight Transform File V1.0\n"
        "#Transform 0\n"
        "Transform: MatrixOffsetTransformBase_double_3_3\n"
        "Parameters: 1 0 0 0 1 0 0 0 1 0 0 0\n"
        "FixedParameters: 0 0 0\n"
    )
