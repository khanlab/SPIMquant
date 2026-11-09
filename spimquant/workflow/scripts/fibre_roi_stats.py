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


def region_table(moments, labels, label_df, metrics):
    """Per-region metrics as a DataFrame (one row per ``label_df`` row).

    Parameters
    ----------
    moments : ndarray (n_moments, z, y, x)
    labels : ndarray int (z, y, x)
    label_df : DataFrame with ``index`` and ``name`` columns
    metrics : list of metric names (keys of ``moments_to_metrics``)
    """
    from fibre_orientation import N_FIBRE, N_TISSUE, aggregate_moments
    from fibre_orientation import moments_to_metrics

    if labels.shape != moments.shape[1:]:
        raise ValueError(
            f"label grid {labels.shape} does not match the moments grid "
            f"{moments.shape[1:]}"
        )
    idx = label_df["index"].to_numpy(dtype=np.int64)
    n_labels = int(max(idx.max(initial=0), labels.max(initial=0))) + 1
    pooled = aggregate_moments(moments, np.clip(labels, 0, None), n_labels)
    values = moments_to_metrics(pooled)
    out = label_df[["index", "name"]].copy()
    for m in metrics:
        out[m] = values[m][idx]
    out["fibrevoxels"] = pooled[N_FIBRE][idx]
    out["tissuevoxels"] = pooled[N_TISSUE][idx]
    return out


def main():
    import dask
    from zarrnii import ZarrNii

    with dask.config.set(scheduler="threads", num_workers=snakemake.threads):
        moments = np.asarray(
            ZarrNii.from_file(snakemake.input.moments).data.compute(),
            dtype=np.float64,
        )
        dseg = ZarrNii.from_nifti(snakemake.input.dseg, axes_order="ZYX")
        labels = np.asarray(dseg.data[0].compute()).astype(np.int64)
    label_df = pd.read_csv(snakemake.input.label_tsv, sep="\t")
    df = region_table(moments, labels, label_df, list(snakemake.params.metrics))
    df.to_csv(snakemake.output.tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
