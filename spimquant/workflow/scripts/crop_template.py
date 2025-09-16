"""
Script to crop template NIfTI files along the X-axis to retain left or right hemisphere.

This script takes a template NIfTI file and crops it along the X-axis (first spatial dimension)
to retain only the left or right hemisphere for registration purposes.
"""

import nibabel as nib
import numpy as np


def crop_template_hemisphere(input_path, output_path, hemisphere):
    """
    Crop a template NIfTI file to retain only left or right hemisphere.

    Parameters:
    -----------
    input_path : str
        Path to input template NIfTI file
    output_path : str
        Path to output cropped NIfTI file
    hemisphere : str
        'left' or 'right' - hemisphere to retain
    """
    # Load the template image
    img = nib.load(input_path)
    data = img.get_fdata()

    # Get the middle point along X-axis (assuming X is the first dimension)
    x_size = data.shape[0]
    x_middle = x_size // 2

    # Crop based on hemisphere
    if hemisphere == "left":
        # Keep the left half (first half of X dimension)
        cropped_data = data[:x_middle, :, :]
    elif hemisphere == "right":
        # Keep the right half (second half of X dimension)
        cropped_data = data[x_middle:, :, :]
    else:
        raise ValueError(f"Hemisphere must be 'left' or 'right', got '{hemisphere}'")

    # Update the affine matrix to reflect the cropping
    affine = img.affine.copy()
    if hemisphere == "right":
        # For right hemisphere, we need to adjust the origin in the affine matrix
        # since we're removing the left part of the image
        affine[0, 3] += x_middle * affine[0, 0]  # Adjust X translation

    # Create new image with cropped data
    cropped_img = nib.Nifti1Image(cropped_data, affine, img.header)

    # Save the cropped image
    nib.save(cropped_img, output_path)


if __name__ == "__main__":
    # Access snakemake variables
    input_template = snakemake.input.template
    output_cropped = snakemake.output.cropped
    hemisphere = snakemake.params.hemisphere

    crop_template_hemisphere(input_template, output_cropped, hemisphere)
