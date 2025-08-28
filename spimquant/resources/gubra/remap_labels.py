import nibabel as nib
import pandas as pd
import json
import numpy as np

# Paths to input files
nifti_path = "../LSFM-mouse-brain-atlas/LSFM_atlas_files/gubra_ano_olf_spacing.nii.gz"  # Path to the input NIfTI file
tsv_path = "../LSFM-mouse-brain-atlas/LSFM_atlas_files/ARA2_annotation_info.csv"  # Path to the TSV file with intensities and names
json_path = "../ABAv3/labelmapper_ABAv3_to_all.json"  # Path to the JSON file with target mappings
output_path = "gubra_ano_olf_spacing_remap.nii.gz"  # Path to save the output NIfTI file

# Load the NIfTI image
nifti_img = nib.load(nifti_path)
nifti_data = nifti_img.get_fdata()

# Load the TSV file into a DataFrame
tsv_data = pd.read_csv(tsv_path)
id_to_name = dict(zip(tsv_data["id"], tsv_data["name"]))

# Load the target mapping from the JSON file
with open(json_path, "r") as f:
    target_map = json.load(f)

# Create a mapping from source intensities to target intensities
intensity_mapping = {}
for entry in target_map:
    output_intensity = entry[0]
    target_name = entry[3]

    # Find the source ID in the TSV file that corresponds to the target name in the JSON
    source_id = next((k for k, v in id_to_name.items() if v == target_name), None)
    if source_id is not None:
        intensity_mapping[source_id] = output_intensity
        print(f"mapping {source_id} to {output_intensity} for {target_name}")
    else:
        print(f"Warning: No match found in TSV for target name '{target_name}' in JSON")

# Remap the intensities in the NIfTI data
remapped_data = np.copy(nifti_data)
for src_intensity, tgt_intensity in intensity_mapping.items():
    remapped_data[nifti_data == src_intensity] = tgt_intensity

# Save the remapped data as a new NIfTI image
remapped_img = nib.Nifti1Image(remapped_data, nifti_img.affine, nifti_img.header)
remapped_img.to_filename(output_path)

print(f"Remapped NIfTI image saved as '{output_path}'")
