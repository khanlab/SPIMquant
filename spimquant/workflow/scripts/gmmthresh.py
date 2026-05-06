"""GMM tail-thresholding segmentation script.

Fits a Gaussian Mixture Model (GMM) with n components to the log-intensity
histogram of a bias-field corrected OME-Zarr image, sorts the components by
mean intensity ascending, then computes a segmentation threshold from the
nth (highest-intensity) component as:

    threshold = expm1(mean_n + k * sigma_n)

where mean_n and sigma_n are derived from the GMM fit in log1p space.
The resulting scalar threshold is applied to the original (linear) image
intensities.  A QC figure showing the histogram, GMM component fits, and
the threshold lines is saved alongside the mask.
"""

import re

import matplotlib
import numpy as np
from dask_setup import get_dask_client
from sklearn.mixture import GaussianMixture
from zarrnii import ZarrNii

matplotlib.use("agg")
import matplotlib.pyplot as plt


def _parse_gmm_method(method: str):
    """Parse a gmm+n{n}k{k} method string.

    k may encode a float using 'p' as the decimal separator, e.g.
    ``n2k2p5`` → n=2, k=2.5.

    Parameters
    ----------
    method:
        Method string such as ``gmm+n2k2p5`` or ``gmm+n3k3``.

    Returns
    -------
    n_components : int
        Number of GMM components.
    k_sigma : float
        Multiplier for σ to set the threshold above the component mean.
    """
    match = re.fullmatch(r"gmm\+n(\d+)k(\d+(?:p\d+)?)", method)
    if match is None:
        raise ValueError(
            f"Invalid gmm seg_method '{method}'. "
            "Expected format: gmm+n<n>k<k> or gmm+n<n>k<int>p<frac>, "
            "e.g. gmm+n2k2 or gmm+n2k2p5"
        )
    n_components = int(match.group(1))
    k_str = match.group(2).replace("p", ".")
    k_sigma = float(k_str)
    return n_components, k_sigma


def _fit_gmm_and_threshold(data_log, n_components, k_sigma, random_state=42):
    """Fit a GMM in log-intensity space and return the threshold and fit info.

    Parameters
    ----------
    data_log:
        1-D array of log1p-transformed intensity values (float32).
    n_components:
        Number of GMM components to fit.
    k_sigma:
        Multiplier applied to the std-dev of the highest-intensity component.
    random_state:
        Random seed for reproducibility.

    Returns
    -------
    threshold_log : float
        Threshold in log-intensity space.
    means_sorted : np.ndarray, shape (n_components,)
        Component means sorted ascending by mean.
    stds_sorted : np.ndarray, shape (n_components,)
        Component standard deviations sorted ascending by mean.
    weights_sorted : np.ndarray, shape (n_components,)
        Component weights sorted ascending by mean.
    """
    if n_components < 1:
        raise ValueError(f"n_components must be >= 1, got {n_components}")

    gmm = GaussianMixture(
        n_components=n_components,
        covariance_type="full",
        random_state=random_state,
        max_iter=200,
    )
    gmm.fit(data_log.reshape(-1, 1))

    # sort components by mean ascending
    order = np.argsort(gmm.means_.ravel())
    means_sorted = gmm.means_.ravel()[order]
    stds_sorted = np.sqrt(gmm.covariances_.ravel()[order])
    weights_sorted = gmm.weights_[order]

    # highest-intensity component (last after ascending sort)
    mean_n = means_sorted[-1]
    sigma_n = stds_sorted[-1]
    threshold_log = mean_n + k_sigma * sigma_n

    return threshold_log, means_sorted, stds_sorted, weights_sorted


def _make_qc_figure(
    bin_edges,
    hist_counts,
    means_sorted,
    stds_sorted,
    weights_sorted,
    threshold_log,
    k_sigma,
    n_components,
):
    """Create a QC figure showing the histogram, GMM fits and threshold lines.

    All quantities are in log1p-intensity space for the GMM overlays; the
    x-axis is labelled in log1p units.
    """
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = bin_edges[1] - bin_edges[0]

    # normalise histogram to a density (area = 1)
    total = hist_counts.sum() * bin_width
    hist_density = hist_counts / (total if total > 0 else 1.0)

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.bar(
        bin_centers,
        hist_density,
        width=bin_width,
        color="steelblue",
        alpha=0.5,
        label="histogram (log1p)",
    )

    # overlay each GMM component
    x_plot = np.linspace(bin_edges[0], bin_edges[-1], 500)
    colors = plt.cm.tab10(np.linspace(0, 1, n_components))
    for idx in range(n_components):
        mu = means_sorted[idx]
        sigma = stds_sorted[idx]
        w = weights_sorted[idx]
        component_pdf = (
            w
            / (sigma * np.sqrt(2 * np.pi))
            * np.exp(-0.5 * ((x_plot - mu) / sigma) ** 2)
        )
        ax.plot(
            x_plot,
            component_pdf,
            color=colors[idx],
            linewidth=1.5,
            label=f"GMM comp {idx + 1}: μ={mu:.3f}, σ={sigma:.3f}, w={w:.3f}",
        )

    # highest-intensity component sigma lines
    mean_n = means_sorted[-1]
    sigma_n = stds_sorted[-1]
    ax.axvline(
        mean_n + sigma_n,
        color="orange",
        linestyle="--",
        linewidth=1.2,
        label=f"μ + σ = {mean_n + sigma_n:.3f}",
    )
    ax.axvline(
        threshold_log,
        color="red",
        linestyle="-",
        linewidth=1.8,
        label=f"threshold (μ + {k_sigma}σ) = {threshold_log:.3f}",
    )

    ax.set_xlabel("log1p(intensity)")
    ax.set_ylabel("density")
    ax.set_title(f"GMM threshold: n={n_components} components, k={k_sigma}")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    with get_dask_client(snakemake.config["dask_scheduler"], snakemake.threads):

        zarrnii_kwargs = snakemake.params.zarrnii_kwargs
        pct_lo, pct_hi = snakemake.params.hist_percentile_range
        bin_width = snakemake.params.hist_bin_width
        n_components = snakemake.params.gmm_n
        k_sigma = snakemake.params.gmm_k

        # ------------------------------------------------------------------ #
        # 1. Load a downsampled version to estimate the percentile-based range
        # ------------------------------------------------------------------ #
        print("estimating intensity range from downsampled image...")
        znimg_ds = None
        for ds_level in [5, 4, 3, 2, 1]:
            try:
                candidate = ZarrNii.from_ome_zarr(
                    snakemake.input.corrected, level=ds_level, **zarrnii_kwargs
                )
                znimg_ds = candidate
                break
            except Exception:
                pass

        if znimg_ds is None:
            znimg_ds = ZarrNii.from_ome_zarr(
                snakemake.input.corrected, **zarrnii_kwargs
            )

        data_ds = znimg_ds.data.compute().ravel().astype(np.float32)

        # validate that we have positive intensities for log transform
        pos_mask = data_ds > 0
        if not pos_mask.any():
            raise ValueError(
                "No positive-intensity voxels found in the downsampled image. "
                "GMM thresholding requires positive intensities."
            )

        range_lo = float(np.percentile(data_ds, pct_lo))
        range_hi = float(np.percentile(data_ds, pct_hi))
        print(
            f"  📊 percentile range [{pct_lo}%, {pct_hi}%]: "
            f"[{range_lo:.3f}, {range_hi:.3f}]"
        )

        # ------------------------------------------------------------------ #
        # 2. Load full-resolution (level=0) corrected image
        # ------------------------------------------------------------------ #
        znimg = ZarrNii.from_ome_zarr(snakemake.input.corrected, **zarrnii_kwargs)

        # ------------------------------------------------------------------ #
        # 3. Compute histogram in linear space using percentile range
        # ------------------------------------------------------------------ #
        n_bins = max(2, int(np.ceil((range_hi - range_lo) / bin_width)))
        print(f"  📊 bins: {n_bins} (bin width: {bin_width})")

        hist_counts, bin_edges = znimg.compute_histogram(
            bins=n_bins, range=[range_lo, range_hi]
        )
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # ------------------------------------------------------------------ #
        # 4. Build a log1p-intensity sample from the histogram for GMM fitting
        # ------------------------------------------------------------------ #
        bin_centers_log = np.log1p(np.maximum(bin_centers, 0.0)).astype(np.float32)

        # expand histogram bins into a sample weighted by counts
        counts_int = hist_counts.astype(np.int64)
        data_log = np.repeat(bin_centers_log, counts_int)

        if data_log.size == 0:
            raise ValueError(
                "Histogram is empty after applying percentile range filter. "
                "Check seg_hist_percentile_range settings."
            )

        print(
            f"fitting GMM with n={n_components} components in log1p space "
            f"(k={k_sigma})..."
        )
        threshold_log, means_sorted, stds_sorted, weights_sorted = (
            _fit_gmm_and_threshold(data_log, n_components, k_sigma)
        )

        # transform threshold back to linear intensity space
        threshold_linear = float(np.expm1(threshold_log))
        print(
            f"  📈 GMM threshold (log1p space): {threshold_log:.4f} "
            f"→ linear: {threshold_linear:.4f}"
        )

        # ------------------------------------------------------------------ #
        # 5. QC figure
        # ------------------------------------------------------------------ #
        bin_edges_log = np.log1p(np.maximum(bin_edges, 0.0))
        hist_counts_log, _ = np.histogram(data_log, bins=bin_edges_log)

        fig = _make_qc_figure(
            bin_edges_log,
            hist_counts_log,
            means_sorted,
            stds_sorted,
            weights_sorted,
            threshold_log,
            k_sigma,
            n_components,
        )
        fig.savefig(snakemake.output.thresholds_png)
        plt.close(fig)

        # ------------------------------------------------------------------ #
        # 6. Apply threshold to original (linear) image and save
        # ------------------------------------------------------------------ #
        print("thresholding image, saving as ome zarr")
        znimg_mask = znimg.segment_threshold(threshold_linear)

        # multiplying binary mask by 100 (so values are 0 and 100) to enable
        # field fraction calculation by subsequent local-mean downsampling
        znimg_mask = znimg_mask * 100

        znimg_mask.to_ome_zarr(snakemake.output.mask, max_layer=5)
