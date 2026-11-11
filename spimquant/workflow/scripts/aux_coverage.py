"""Fraction of each atlas region covered by an additional scan's field of view.

The field of view is the additional scan's nonzero voxels (zero-padded margins
hold no data). It is resampled (linearly, so partial voxels count partly) into
the reference scan's atlas labels, which cover the whole sample, and averaged
over each region. Only the index, name and coverage columns are written, so the
table merges into mergedsegstats.tsv without clashing with its columns.
"""

import os
import subprocess
import tempfile

import nibabel as nib
import numpy as np
from zarrnii import ZarrNii, ZarrNiiAtlas


def main(snakemake):
    fov_img = nib.load(snakemake.input.fov)
    fov = (np.asarray(fov_img.dataobj) > 0).astype(np.float32)

    with tempfile.TemporaryDirectory() as tmp:
        fov_path = os.path.join(tmp, "fov.nii.gz")
        fov_in_ref = os.path.join(tmp, "fov_in_ref.nii.gz")
        nib.save(nib.Nifti1Image(fov, fov_img.affine), fov_path)
        subprocess.run(
            [
                "greedy",
                "-threads",
                str(snakemake.threads),
                "-d",
                "3",
                "-rf",
                snakemake.input.dseg,
                "-rm",
                fov_path,
                fov_in_ref,
                "-r",
                *snakemake.params.to_aux.split(),
            ],
            check=True,
        )
        atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)
        df = atlas.aggregate_image_by_regions(
            ZarrNii.from_file(fov_in_ref),
            aggregation_func="mean",
            column_name="coverage",
        )

    df = df[["index", "name", "coverage"]]
    df["coverage"] = df["coverage"].clip(0, 1)
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)


if __name__ == "__main__":
    main(snakemake)  # noqa: F821
