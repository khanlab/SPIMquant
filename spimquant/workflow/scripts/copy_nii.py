import shutil

src = str(snakemake.input[0])
dst = str(snakemake.output[0])


def _ext(path):
    return ".nii.gz" if path.endswith(".nii.gz") else ".nii"


if _ext(src) == _ext(dst):
    # same format: byte copy, avoids decompressing the whole volume into
    # memory (localrules run on the submit host, where large templates
    # can exceed the login node's memory limits and get SIGKILLed)
    shutil.copyfile(src, dst)
else:
    import nibabel as nib

    nib.load(src).to_filename(dst)
