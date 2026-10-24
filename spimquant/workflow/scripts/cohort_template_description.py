"""Write template_description.json and a SPIMquant templates: config snippet.

The snippet points at Zenodo record URLs when a deposition ID was given, and
otherwise at paths relative to spimquant/resources/ (after copying the
tpl-{name} directory there).
"""

import json
import os
from datetime import date

import nibabel as nib
import yaml

p = snakemake.params
name = p.name
orig = p.orig_template_config

files = {os.path.basename(f): f for f in snakemake.input.files}


def ref(fname):
    if fname not in files:
        raise ValueError(f"{fname} is not one of the exported files")
    if p.zenodo_deposition is not None:
        return f"https://zenodo.org/records/{p.zenodo_deposition}/files/{fname}"
    return f"tpl-{name}/{fname}"


anat_fname = f"tpl-{name}_anat.nii.gz"
entry = {"anat": ref(anat_fname)}
if f"tpl-{name}_desc-brain_mask.nii.gz" in files:
    entry["mask"] = ref(f"tpl-{name}_desc-brain_mask.nii.gz")
# stain order is the registration priority SPIMquant uses
entry["spim_templates"] = {
    stain: ref(f"tpl-{name}_stain-{stain}_SPIM.nii.gz") for stain in p.stains
}
if "default_segs" in orig:
    entry["default_segs"] = list(orig["default_segs"])
entry["atlases"] = {
    seg: {
        "dseg": ref(f"tpl-{name}_seg-{seg}_dseg.nii.gz"),
        "tsv": ref(f"tpl-{name}_seg-{seg}_dseg.tsv"),
    }
    for seg in orig.get("atlases", {})
}

if p.zenodo_deposition is not None:
    location = (
        f"# Files are in Zenodo deposition {p.zenodo_deposition}; the URLs resolve\n"
        "# once the deposition is published.\n"
    )
else:
    location = (
        f"# Paths are relative to spimquant/resources/: copy {p.template_dir}\n"
        "# there, or rerun cohort_template_export with\n"
        "# --cohort-template-zenodo-deposition to upload to Zenodo.\n"
    )

with open(snakemake.output.config_snippet, "w") as f:
    f.write(
        f"# SPIMquant templates: entry for tpl-{name}, a cohort template built from\n"
        f"# {len(p.subjects)} subjects registered to tpl-{p.orig_template}, with atlases\n"
        f"# warped from tpl-{p.orig_template}.\n"
        f"# Add it under templates: in spimquant/config/snakebids.yml, and add {name}\n"
        "# to the --template choices.\n" + location
    )
    yaml.safe_dump({"templates": {name: entry}}, f, sort_keys=False)


anat = nib.load(files[anat_fname])
description = {
    "Name": f"SPIMquant cohort template {name}",
    "Identifier": name,
    "Species": "Mus musculus",
    "DateCreated": date.today().isoformat(),
    "GeneratedBy": [
        {
            "Name": "SPIMquant",
            "Description": "cohort_template analysis level: FireANTs multivariate "
            "template building, then greedy registration of the original template "
            "to the cohort template to warp its atlases",
            "CodeURL": "https://github.com/khanlab/SPIMquant",
        },
        {
            "Name": "FireANTs",
            "CodeURL": "https://github.com/rohitrango/FireANTs",
        },
    ],
    "OriginalTemplate": p.orig_template,
    "Stains": list(p.stains),
    "InputLevel": p.level,
    # isotropic voxel size (mm) the inputs were resampled to, null if not resampled
    "InputResampleMM": p.resample_mm,
    "TemplateIterations": p.iterations,
    "Subjects": list(p.subjects),
    "Atlases": list(entry["atlases"].keys()),
    "Transforms": {
        f"tpl-{name}_from-{p.orig_template}_to-{name}_xfm.nii.gz": f"resamples "
        f"tpl-{p.orig_template} images into tpl-{name} space",
        f"tpl-{name}_from-{name}_to-{p.orig_template}_xfm.nii.gz": f"resamples "
        f"tpl-{name} images into tpl-{p.orig_template} space",
    },
    "FireANTsParameters": dict(p.fireants),
    "RegistrationParameters": dict(p.reg),
    "res": {
        "1": {
            "shape": list(anat.shape[:3]),
            "zooms": [float(z) for z in anat.header.get_zooms()[:3]],
            "origin": [float(o) for o in anat.affine[:3, 3]],
        }
    },
}

with open(snakemake.output.description, "w") as f:
    json.dump(description, f, indent=2)
