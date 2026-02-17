"""Average multiple registered MRI images with optional upsampling.

This script takes multiple MRI images that have been N4-corrected and registered
to a common reference, applies their transformations with optional upsampling,
and averages them to create a super-resolved or higher SNR final image.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np


def apply_transform_with_upsample(
    input_nii, xfm_ras, ref_nii, output_nii, upsample_factor=1, threads=1
):
    """Apply transformation to an image with optional upsampling.

    Args:
        input_nii: Path to input image
        xfm_ras: Path to transformation matrix
        ref_nii: Path to reference image
        output_nii: Path to output image
        upsample_factor: Factor to upsample the reference grid (default: 1)
        threads: Number of threads for greedy
    """
    if upsample_factor == 1:
        # No upsampling, direct application
        cmd = [
            "greedy",
            "-threads",
            str(threads),
            "-d",
            "3",
            "-rf",
            ref_nii,
            "-rm",
            input_nii,
            output_nii,
            "-r",
            xfm_ras,
        ]
        subprocess.run(cmd, check=True)
    else:
        # Create upsampled reference grid
        ref_img = nib.load(ref_nii)
        new_shape = tuple(int(s * upsample_factor) for s in ref_img.shape)

        # Create new affine with scaled voxel sizes
        # Get original voxel sizes (magnitude of each affine column)
        orig_affine = ref_img.affine
        new_affine = orig_affine.copy()

        # Scale each spatial dimension by dividing voxel size
        # This is done by dividing each column vector by the upsample factor
        for i in range(3):
            new_affine[:3, i] = orig_affine[:3, i] / upsample_factor

        # Create temporary upsampled reference
        with tempfile.NamedTemporaryFile(suffix=".nii.gz", delete=False) as tmp_ref:
            upsampled_ref = nib.Nifti1Image(
                np.zeros(new_shape, dtype=np.float32), new_affine
            )
            nib.save(upsampled_ref, tmp_ref.name)
            tmp_ref_path = tmp_ref.name

        try:
            # Apply transform with upsampled reference
            cmd = [
                "greedy",
                "-threads",
                str(threads),
                "-d",
                "3",
                "-rf",
                tmp_ref_path,
                "-rm",
                input_nii,
                output_nii,
                "-r",
                xfm_ras,
            ]
            subprocess.run(cmd, check=True)
        finally:
            Path(tmp_ref_path).unlink(missing_ok=True)


def main():
    # Access snakemake inputs
    ref_nii = snakemake.input.ref
    mri_files = snakemake.input.mri
    xfm_files = snakemake.input.xfm
    output_nii = snakemake.output.nii
    upsample_factor = snakemake.params.upsample_factor
    threads = snakemake.threads

    # Handle single MRI case
    if len(mri_files) == 1:
        # Just copy the single MRI to output (it's already N4 corrected)
        shutil.copy(mri_files[0], output_nii)
        return

    # Transform all MRIs to reference space with optional upsampling
    transformed_files = []
    temp_dir = Path(tempfile.mkdtemp())

    try:
        for idx, (mri, xfm) in enumerate(zip(mri_files, xfm_files)):
            temp_output = temp_dir / f"transformed_{idx}.nii.gz"
            apply_transform_with_upsample(
                mri, xfm, ref_nii, str(temp_output), upsample_factor, threads
            )
            transformed_files.append(str(temp_output))

        # Load all transformed images
        imgs = [nib.load(f) for f in transformed_files]
        data_arrays = [img.get_fdata() for img in imgs]

        # Average the images
        averaged_data = np.mean(data_arrays, axis=0)

        # Save the averaged image using the first image's header/affine
        averaged_img = nib.Nifti1Image(averaged_data.astype(np.float32), imgs[0].affine)
        nib.save(averaged_img, output_nii)

    finally:
        # Clean up temporary files
        for f in transformed_files:
            Path(f).unlink(missing_ok=True)
        try:
            temp_dir.rmdir()
        except OSError:
            # Directory may not be empty, use shutil
            shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
