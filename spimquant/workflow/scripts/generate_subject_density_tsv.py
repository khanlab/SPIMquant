import nibabel as nib
import pandas as pd
import numpy as np

dseg_df = pd.read_csv(snakemake.input.volumes_tsv,sep='\t')

blobs_df = pd.read_csv(snakemake.input.blobs_tsv,sep='\t')

merged_df = dseg_df.merge(blobs_df.groupby('label_index').size().reset_index(name='blob_count'), how='left', left_on='index', right_on='label_index')

merged_df = merged_df.drop(columns=['label_index'])

merged_df['density'] = merged_df['blob_count'] / merged_df['volume']

merged_df.to_csv(snakemake.output.density_tsv, sep="\t", index=False)
