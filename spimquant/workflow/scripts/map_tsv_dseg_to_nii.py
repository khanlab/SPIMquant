from zarrnii import ZarrNiiAtlas
import pandas as pd

atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

feature_data = pd.read_csv(snakemake.input.tsv, sep="\t")

img = atlas.create_feature_map(
    feature_data,
    feature_column=snakemake.params.feature_column,
    label_column=snakemake.params.label_column,
)

#force atlas units (until zarrnii issue #203 fixed):
img.ngff_image.axes_units = {'x': 'millimeter', 'y': 'millimeter', 'z': 'millimeter'}

img.to_nifti(snakemake.output.nii)
