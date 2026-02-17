"""Copy preprocessed MRI or apply N4 to single MRI.

For multi-MRI cases, the input is already N4-corrected and averaged (desc-preproc),
so we just copy it. For single MRI cases, we apply N4 correction.
"""

import shutil
import subprocess

if __name__ == "__main__":
    input_nii = snakemake.input.nii
    output_nii = snakemake.output.nii
    is_multi_mri = snakemake.params.is_multi_mri

    if is_multi_mri:
        # Multi-MRI case: input is already N4-corrected, just copy
        shutil.copy(input_nii, output_nii)
    else:
        # Single MRI case: apply N4
        cmd = [
            "N4BiasFieldCorrection",
            "-i",
            input_nii,
            "-o",
            output_nii,
            "-d",
            "3",
            "-v",
        ]
        subprocess.run(cmd, check=True)
