"""Distances between the objects of all stains (rule annotate_object_distances).

For every stain S with objects, each object (of any stain) gets:

* ``S_dist``: distance from its centroid to the nearest centroid of an object of
  S. For objects of S itself, the nearest *other* object (nearest-neighbour
  spacing).
* ``S_edgedist``: ``S_dist`` minus the equivalent-sphere radius of that nearest
  object (from its ``volume``), clipped at 0. An approximate distance to the
  nearest object's surface: an object touching a large plaque reads ~0 rather
  than the plaque's radius.

Distances use the subject-space physical centroids (mm). Sampled signed
distance columns (``<stain>_sdt``, one parquet each, row-aligned) are joined
on.
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def equivalent_radius(volume):
    """Radius of the sphere of the given volume."""
    return np.cbrt(3.0 * np.asarray(volume, dtype=np.float64) / (4.0 * np.pi))


def add_object_distances(df, coord_cols, stain_col="stain", volume_col="volume"):
    """Return ``df`` with ``<stain>_dist`` and ``<stain>_edgedist`` columns.

    Stains with a single object have no "other" object for their own rows,
    so those rows get NaN; a stain's rows are NaN against an empty stain.
    """
    df = df.copy()
    coords = df[coord_cols].to_numpy(dtype=np.float64)
    labels = df[stain_col].to_numpy()
    radius = (
        equivalent_radius(df[volume_col])
        if volume_col in df.columns
        else np.zeros(len(df))
    )
    for stain in pd.unique(labels):
        sel = np.flatnonzero(labels == stain)
        dist = np.full(len(df), np.nan)
        edge = np.full(len(df), np.nan)
        if len(sel) > 0:
            tree = cKDTree(coords[sel])
            own = labels == stain
            other = ~own
            if other.any():
                d, i = tree.query(coords[other], k=1)
                dist[other] = d
                edge[other] = np.maximum(d - radius[sel[i]], 0.0)
            if own.sum() > 1:
                # k=2: the nearest hit is the object itself
                d, i = tree.query(coords[own], k=2)
                dist[own] = d[:, 1]
                edge[own] = np.maximum(d[:, 1] - radius[sel[i[:, 1]]], 0.0)
        df[f"{stain}_dist"] = dist
        df[f"{stain}_edgedist"] = edge
    return df


def main():
    df = pd.read_parquet(snakemake.input.regionprops_parquet)
    coord_cols = list(snakemake.params.coord_column_names)
    if len(df) > 0:
        df = add_object_distances(df, coord_cols)
    for path in snakemake.input.sdt_parquets:
        sdt = pd.read_parquet(path)
        if len(sdt) != len(df):
            raise ValueError(
                f"{path} has {len(sdt)} rows but the objects table has {len(df)}"
            )
        for col in sdt.columns:
            df[col] = sdt[col].to_numpy()
    df.to_parquet(snakemake.output.regionprops_parquet, index=False)


if __name__ == "__main__" and "snakemake" in dir():
    main()
