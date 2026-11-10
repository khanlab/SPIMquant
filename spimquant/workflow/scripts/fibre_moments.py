"""Fibre orientation moments of a bias-corrected SPIM image (rule fibre_moments).

Runs the multiscale Frangi filter chunk by chunk and reduces each chunk to
per-super-voxel moments (see ``fibre_orientation.py``), written as an
OME-Zarr with one channel per moment.  The Frangi background parameter
``gamma`` and the fibre threshold are estimated first from voxels sampled in
tissue-rich chunks (chosen on a coarser pyramid level, as in dogthresh), and
saved with the other parameters to a JSON file, alongside a QC histogram.
"""

import json

import numpy as np

_RANDOM_STATE = 42


def block_factor(block_um, spacing_um, chunk=None):
    """Super-voxel size in voxels per axis for an edge of ``block_um``.

    With ``chunk``, each factor is snapped (in log scale) to the nearest
    divisor of ``chunk``, so that chunks stay aligned with the stored ones
    and the input needs no rechunking (~12% faster on a mouse brain).
    """
    out = []
    for s in spacing_um:
        f = max(1.0, float(block_um) / float(s))
        if chunk:
            divisors = [d for d in range(1, int(chunk) + 1) if int(chunk) % d == 0]
            f = min(divisors, key=lambda d: abs(np.log(d / f)))
        out.append(max(1, int(round(f))))
    return tuple(out)


def load_lowres(path, level_offset):
    """Lowres pyramid level ``level_offset`` of ``path`` (or the coarsest)."""
    import dask
    from zarrnii import ZarrNii

    with dask.config.set(scheduler="threads"):
        level = int(level_offset)
        while True:
            try:
                znimg = ZarrNii.from_file(path, level=level)
                break
            except Exception:
                if level == 0:
                    raise
                level -= 1
        return np.asarray(znimg.data[0].compute()), level


def choose_blocks(lowres, chunks, shape, n_blocks):
    """Tissue-rich chunk indices (as dogthresh chooses its stat blocks)."""
    from dogthresh import select_stat_blocks, tissue_fraction_per_block

    fractions = tissue_fraction_per_block(lowres, chunks, shape)
    blocks = select_stat_blocks(fractions, n_blocks, seed=_RANDOM_STATE)
    if not blocks:
        raise ValueError("No tissue-containing chunks found for Frangi statistics")
    return blocks


def make_qc_figure(est, scales_um, method):
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    from fibre_orientation import structureness

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    vals = est["frangi_samples"]
    ax.hist(vals, bins=200, range=(0, 1), color="steelblue", alpha=0.7)
    ax.set_yscale("log")
    thr = est["threshold"]
    frac = float(np.mean(vals > thr)) if vals.size else 0.0
    ax.axvline(
        thr,
        color="red",
        linewidth=1.8,
        label=f"threshold ({method}) = {thr:.3f}\n{100 * frac:.2f}% of sampled tissue",
    )
    ax.set_xlabel("multiscale Frangi response")
    ax.set_ylabel("sampled tissue voxels (log)")
    ax.legend(fontsize=8)
    ax.set_title("Fibre threshold")

    ax = axes[1]
    for s, (scale, w, g) in enumerate(
        zip(scales_um, est["eig_samples"], est["gammas"])
    ):
        sn = structureness(w)
        sn = sn[np.isfinite(sn)]
        if sn.size == 0:
            continue
        hi = np.percentile(sn, 99.99)
        (line,) = ax.plot([], [])
        ax.hist(
            sn,
            bins=200,
            range=(0, hi),
            histtype="step",
            color=line.get_color(),
            label=f"sigma {scale:g} um (gamma {g:.3g})",
        )
        ax.axvline(g, color=line.get_color(), linestyle="--", linewidth=1)
    ax.set_yscale("log")
    ax.set_xlabel("scale-normalised Hessian norm")
    ax.set_title("Frangi gamma per scale")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def main():
    import tempfile
    from contextlib import ExitStack

    from dask_setup import get_dask_client
    from fibre_orientation import aligned_chunks, n_moments, zarrnii_fibre_moments
    from zarrnii import ZarrNii

    params = snakemake.params
    scales_um = [float(s) for s in params.scales_um]

    with ExitStack() as stack:
        stack.enter_context(
            get_dask_client(
                snakemake.config["dask_scheduler"],
                snakemake.threads,
                threads_per_worker=int(params.dask_threads_per_worker),
                mem_mb=snakemake.resources.get("mem_mb"),
                reserve_mb=float(params.dask_reserve_mb),
            )
        )
        znimg = ZarrNii.from_file(snakemake.input.corrected)
        if znimg.data.ndim != 4:
            raise ValueError(
                f"Expected a (c, z, y, x) image, got shape {znimg.data.shape}"
            )
        spacing_um = [z * 1000.0 for z in znimg.get_zooms()]  # mm -> um
        chunk = int(params.chunk_size)
        factor = block_factor(params.block_um, spacing_um, chunk)
        shape = znimg.data.shape[1:]
        chunks = aligned_chunks(shape, (chunk,) * 3, factor)
        print(
            f"voxel size (um) {np.round(spacing_um, 3)}, super-voxel {factor} vox "
            f"({np.round(np.multiply(factor, spacing_um), 1)} um), "
            f"chunks {[c[0] for c in chunks]}, scales {scales_um} um"
        )

        lowres, lowres_level = load_lowres(
            snakemake.input.corrected, params.stat_level_offset
        )
        block_ids = choose_blocks(lowres, chunks, shape, int(params.stat_blocks))
        del lowres
        print(f"sampling {len(block_ids)} chunks (chosen on level +{lowres_level})")

        moments, est = zarrnii_fibre_moments(
            znimg,
            scales_um,
            factor,
            lmax=int(params.lmax),
            block_ids=block_ids,
            chunk_size=(chunk,) * 3,
            alpha=float(params.alpha),
            beta=float(params.beta),
            edge_min_weight=float(params.edge_min_weight),
            n_samples_per_block=int(params.stat_samples_per_block),
            gamma_quantile=float(params.gamma_quantile),
            gamma_factor=float(params.gamma_factor),
            threshold_method=params.threshold_method,
            seed=_RANDOM_STATE,
            # computed slab by slab into a temporary zarr: bounds the halo
            # data dask holds on the workers (see store_moments_by_slab)
            slab_store=stack.enter_context(
                tempfile.TemporaryDirectory(
                    suffix="_moments.zarr", dir=snakemake.resources.get("tmpdir")
                )
            ),
            slab_chunks=int(params.slab_chunks),
        )
        print("moments computed")
        gammas = [float(g) for g in est["gammas"]]
        print(f"gammas {np.round(gammas, 4)}, threshold {est['threshold']:.4f}")

        fig = make_qc_figure(est, scales_um, params.threshold_method)
        fig.savefig(snakemake.output.qc_png)

        with open(snakemake.output.params_json, "w") as f:
            json.dump(
                {
                    "scales_um": scales_um,
                    "spacing_um": [float(s) for s in spacing_um],
                    "block_factor": list(factor),
                    "block_um": [float(f * s) for f, s in zip(factor, spacing_um)],
                    "lmax": int(params.lmax),
                    "n_moments": n_moments(int(params.lmax)),
                    "alpha": float(params.alpha),
                    "beta": float(params.beta),
                    "gammas": gammas,
                    "gamma_quantile": float(params.gamma_quantile),
                    "gamma_factor": float(params.gamma_factor),
                    "threshold": float(est["threshold"]),
                    "threshold_method": params.threshold_method,
                    "edge_min_weight": float(params.edge_min_weight),
                    "stat_blocks": [list(map(int, b)) for b in block_ids],
                    "sampled_voxels": int(est["frangi_samples"].size),
                },
                f,
                indent=2,
            )

        print("saving moments as ome zarr")
        # moments are sums, so a mean-downsampled pyramid would be misleading
        moments.to_ome_zarr(
            snakemake.output.moments,
            max_layer=0,
            **snakemake.config["zarrnii_out_kwargs"],
        )


if __name__ == "__main__":
    main()
