from zarrnii import ZarrNii, ZarrNiiAtlas
import numpy as np
import pandas as pd

centroids = np.load(snakemake.input.centroids_npy)

atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

df_centroids, df_counts = atlas.label_centroids(centroids, include_names=True)

df_centroids.to_csv(snakemake.output.centroids_tsv, sep='\t', index=False)
df_counts.to_csv(snakemake.output.counts_tsv, sep='\t', index=False)
