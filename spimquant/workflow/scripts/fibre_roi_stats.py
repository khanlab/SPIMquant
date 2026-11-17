"""Atlas-region fibre orientation statistics from pooled moments.

The fibre moments (counts, orientation tensor and SH sums per super-voxel)
are summed over each atlas region (labels resampled onto the super-voxel
grid), and the metrics are computed from the pooled sums.  This is the
region's own dispersion (e.g. of all its fibres' orientations), not an
average of local dispersion values.

Writes one row per label in the label TSV: ``index``, ``name``, the metrics
in ``params.metrics``, and the pooled ``fibrevoxels`` / ``tissuevoxels``.
"""

import numpy as np
import pandas as pd


def pooled_moments(moments, labels, n_labels, slab_bytes=5e8):
    """Moments summed per label, reading (c, z, y, x) ``moments`` in z-slabs.

    ``moments`` and ``labels`` may be dask or numpy arrays on the same grid.
    Returns float64 (n_moments, n_labels).
    """
    from fibre_orientation import aggregate_moments

    if tuple(labels.shape) != tuple(moments.shape[1:]):
        raise ValueError(
            f"label grid {tuple(labels.shape)} does not match the moments grid "
            f"{tuple(moments.shape[1:])}"
        )
    c, nz, ny, nx = moments.shape
    step = max(1, int(slab_bytes // (8 * c * ny * nx)))
    pooled = np.zeros((c, n_labels), dtype=np.float64)
    for z0 in range(0, nz, step):
        z1 = min(z0 + step, nz)
        m = np.asarray(moments[:, z0:z1], dtype=np.float64)
        lab = np.clip(np.asarray(labels[z0:z1]).astype(np.int64), 0, None)
        pooled += aggregate_moments(m, lab, n_labels)
    return pooled


def region_table(moments, labels, label_df, metrics, max_label=None, min_fibre=0.0):
    """Per-region metrics as a DataFrame (one row per ``label_df`` row).

    Parameters
    ----------
    moments : array (n_moments, z, y, x), numpy or dask
    labels : array int (z, y, x), numpy or dask
    label_df : DataFrame with ``index`` and ``name`` columns
    metrics : list of metric names (keys of ``moments_to_metrics``)
    max_label : int, optional
        Largest label value in ``labels`` (computed if not given).
    min_fibre : float
        Orientation metrics are NaN for regions with at most this many
        (weighted) fibre voxels.
    """
    from fibre_orientation import N_FIBRE, N_TISSUE, moments_to_metrics

    idx = label_df["index"].to_numpy(dtype=np.int64)
    if max_label is None:
        max_label = int(np.asarray(labels).max(initial=0))
    n_labels = int(max(idx.max(initial=0), max_label)) + 1
    pooled = pooled_moments(moments, labels, n_labels)
    values = moments_to_metrics(pooled, min_fibre=min_fibre, missing=np.nan)
    out = label_df[["index", "name"]].copy()
    for m in metrics:
        out[m] = values[m][idx]
    out["fibrevoxels"] = pooled[N_FIBRE][idx]
    out["tissuevoxels"] = pooled[N_TISSUE][idx]
    return out


def main():
    import dask
    from fibre_orientation import min_fibre_from_params
    from zarrnii import ZarrNii

    min_fibre = min_fibre_from_params(
        snakemake.input.params_json, snakemake.params.min_fraction
    )
    with dask.config.set(scheduler="threads", num_workers=snakemake.threads):
        moments = ZarrNii.from_file(snakemake.input.moments).data
        dseg = ZarrNii.from_nifti(snakemake.input.dseg, axes_order="ZYX")
        # labels on the super-voxel grid are small (1/n_moments of moments)
        labels = np.asarray(dseg.data[0].compute()).astype(np.int64)
        label_df = pd.read_csv(snakemake.input.label_tsv, sep="\t")
        df = region_table(
            moments,
            labels,
            label_df,
            list(snakemake.params.metrics),
            min_fibre=min_fibre,
        )
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
