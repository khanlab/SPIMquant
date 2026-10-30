"""Write template_description.json and a SPIMquant templates: config snippet.

The snippet is a SPIMquant config file that selects the cohort template, for
use with --configfile. It adds a template source: the Zenodo record when a
deposition ID was given, and otherwise the exported files' directory.
"""

import json
import os
from datetime import date

import nibabel as nib
import yaml

p = snakemake.params
name = p.name
files = {os.path.basename(f): f for f in snakemake.input.files}
anat_fname = f"tpl-{name}_anat.nii.gz"

# the exported files are named so the template components find them
if p.zenodo_deposition is not None:
    source = {"type": "zenodo", "record": p.zenodo_deposition}
else:
    source = {"type": "local", "path": os.path.abspath(p.template_dir)}

# stain order is the registration priority SPIMquant uses
settings = {"spim_stains": list(p.stains)}
if "default_segs" in p.orig_template_settings:
    settings["default_segs"] = list(p.orig_template_settings["default_segs"])

if p.zenodo_deposition is not None:
    location = (
        f"# Files are in Zenodo deposition {p.zenodo_deposition}; the URLs resolve\n"
        "# once the deposition is published.\n"
    )
else:
    location = (
        f"# Paths point at the exported files in {os.path.abspath(p.template_dir)};\n"
        "# rerun cohort_template_export with --cohort-template-zenodo-deposition to\n"
        "# upload to Zenodo instead.\n"
    )

with open(snakemake.output.config_snippet, "w") as f:
    f.write(
        f"# SPIMquant config for tpl-{name}, a cohort template built from\n"
        f"# {len(p.subjects)} subjects registered to tpl-{p.orig_template}, with atlases\n"
        f"# warped from tpl-{p.orig_template}. Use it with:\n"
        f"#   spimquant ... --configfile tpl-{name}_spimquant.yml\n"
        "# To add it to SPIMquant itself, move the source and template_settings\n"
        f"# entries into spimquant/config/snakebids.yml and add {name} to the\n"
        "# --template choices.\n" + location
    )
    yaml.safe_dump(
        {
            "template": name,
            "templates": {"sources": {f"tpl-{name}": source}},
            "template_settings": {name: settings},
        },
        f,
        sort_keys=False,
    )


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
    "Atlases": list(p.atlas_segs),
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
