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

    if dst.endswith(".gz"):
        # verify the gzip stream is complete (streamed, so no full decompression
        # in memory), since a truncated download otherwise reads without error
        # as a partially-zeroed image
        import gzip

        try:
            with gzip.open(dst, "rb") as f:
                while f.read(1 << 24):
                    pass
        except (EOFError, OSError) as e:
            import os

            os.remove(dst)
            raise RuntimeError(
                f"{src} is not a complete gzip file ({e}); it may be a truncated "
                "download, delete it and re-run"
            ) from e
else:
    import nibabel as nib

    nib.load(src).to_filename(dst)
