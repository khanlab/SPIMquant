"""Map computed centroids to atlas regions and generate labeled statistics.

This script takes computed centroids from segmentation and maps them to atlas
regions, generating two outputs:
1. Labeled centroids with their corresponding atlas region assignments
2. Count statistics showing the number of centroids per atlas region
"""
from zarrnii import ZarrNii, ZarrNiiAtlas
import numpy as np
import pandas as pd

centroids = np.load(snakemake.input.centroids_npy)

atlas = ZarrNiiAtlas.from_files(snakemake.input.dseg, snakemake.input.label_tsv)

df_centroids, df_counts = atlas.label_centroids(centroids, include_names=True)

df_centroids.to_csv(snakemake.output.centroids_tsv, sep="\t", index=False)
df_counts.to_csv(snakemake.output.counts_tsv, sep="\t", index=False)
