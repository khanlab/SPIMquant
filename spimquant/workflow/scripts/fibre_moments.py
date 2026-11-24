"""Fibre orientation moments of a bias-corrected SPIM image (rule fibre_moments).

Runs the multiscale Frangi filter chunk by chunk and reduces each chunk to
per-super-voxel moments (see ``fibre_orientation.py``), written as an
OME-Zarr with one channel per moment.  The image is filtered at one pyramid
level and the moments are written at a coarser one: each voxel of the output
level is one super-voxel (``level_factor``).  The Frangi scales follow the
voxel size of the filtered image (``resolve_scales``) unless given in
micrometres.  The Frangi background
parameter ``gamma`` and the fibre threshold are estimated first from voxels
sampled in tissue-rich chunks (chosen on a coarser pyramid level, as in
dogthresh), and saved with the other parameters to a JSON file, alongside a
QC histogram.
"""

import json

import numpy as np

_RANDOM_STATE = 42


def resolve_scales(
    spacing_um, scales_um=None, max_factor=2.0, min_factor=1.0, per_octave=2.0
):
    """Frangi scales (um): ``scales_um`` if given, else from the voxel size
    (``fibre_orientation.frangi_scales``)."""
    from fibre_orientation import frangi_scales

    if scales_um:
        return [float(s) for s in scales_um]
    return frangi_scales(spacing_um, max_factor, min_factor, per_octave)


# a super-voxel should span several of the largest structures detected
_MIN_BLOCK_SCALES = 3.0


def level_factor(path, level_offset):
    """Voxels of ``path`` per voxel of its pyramid level ``level_offset``,
    per spatial axis.

    The super-voxels are the voxels of that coarser level, so the factor is
    the ratio of the two voxel sizes.  Pyramids written by zarrnii keep a
    coarse axis until the others catch up (near-isotropic), so the factor can
    differ between axes: 3.25 x 3.25 x 5.5 um voxels and an offset of 3 give
    8 x 8 x 4 voxels (26 x 26 x 22 um).  A level beyond the stored ones is
    derived by zarrnii from the coarsest stored one.
    """
    from zarrnii import ZarrNii

    level_offset = int(level_offset)
    if level_offset < 0:
        raise ValueError("the fibre output level cannot be finer than the fibre level")
    if level_offset == 0:
        return (1, 1, 1)
    fine = ZarrNii.from_file(path).get_zooms()
    coarse = ZarrNii.from_file(path, level=level_offset).get_zooms()
    return tuple(max(1, int(round(c / f))) for c, f in zip(coarse, fine))


def block_warning(factor, spacing_um, scales_um):
    """A warning if the super-voxels are small for the largest scale, else None."""
    edge = min(f * s for f, s in zip(factor, spacing_um))
    need = _MIN_BLOCK_SCALES * max(scales_um)
    if edge >= need:
        return None
    return (
        f"WARNING: super-voxels of {edge:.3g} um are less than "
        f"{_MIN_BLOCK_SCALES:g}x the largest scale ({max(scales_um):.3g} um): "
        "the local orientation metrics will describe single structures; use a "
        "coarser --fibre_output_level"
    )


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


def threshold_label(method, factor, quantile):
    """How the fibre threshold was set, for the QC figure."""
    if method == "quantile":
        return f"{float(factor):g} x q{100 * float(quantile):g}"
    return str(method)


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
        scales_um = resolve_scales(
            spacing_um,
            params.scales_um,
            params.scale_max_factor,
            params.scale_min_factor,
            params.scales_per_octave,
        )
        filter_level = int(params.filter_level)
        output_level = int(snakemake.wildcards.level)
        chunk = int(params.chunk_size)
        factor = level_factor(snakemake.input.corrected, output_level - filter_level)
        if any(chunk % f for f in factor):
            raise ValueError(
                f"super-voxel {factor} voxels does not divide fibre_chunk_size "
                f"({chunk}); use a finer --fibre_output_level or a larger chunk"
            )
        warning = block_warning(factor, spacing_um, scales_um)
        if warning:
            print(warning, flush=True)
        shape = znimg.data.shape[1:]
        chunks = aligned_chunks(shape, (chunk,) * 3, factor)
        print(
            f"level {filter_level} -> {output_level}: "
            f"voxel size (um) {np.round(spacing_um, 3)}, super-voxel {factor} vox "
            f"({np.round(np.multiply(factor, spacing_um), 1)} um), "
            f"chunks {[c[0] for c in chunks]}, scales {np.round(scales_um, 2)} um"
        )

        lowres, lowres_level = load_lowres(
            snakemake.input.corrected, params.stat_level_offset
        )
        block_ids = choose_blocks(lowres, chunks, shape, int(params.stat_blocks))
        del lowres
        print(f"sampling {len(block_ids)} chunks (chosen on level +{lowres_level})")

        psf = params.psf_fwhm_um
        if isinstance(psf, (list, tuple)):
            psf = [float(v) for v in psf]  # (z, y, x) in um
        kwargs = dict(
            scales_um=scales_um,
            factor=factor,
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
            threshold_factor=float(params.threshold_factor),
            threshold_quantile=float(params.threshold_quantile),
            seed=_RANDOM_STATE,
            # computed slab by slab into a temporary zarr: bounds the halo
            # data dask holds on the workers (see store_moments_by_slab)
            slab_store=stack.enter_context(
                tempfile.TemporaryDirectory(
                    suffix="_moments.zarr", dir=snakemake.resources.get("tmpdir")
                )
            ),
            slab_chunks=int(params.slab_chunks),
            soft_lo=float(params.soft_lo),
            local_norm_um=(
                None if params.local_norm_um is None else float(params.local_norm_um)
            ),
        )
        try:
            moments, est = zarrnii_fibre_moments(znimg, psf_fwhm_um=psf, **kwargs)
        except ValueError as e:
            if psf != "auto":
                raise
            print(f"WARNING: PSF estimate failed ({e}); no PSF equalisation")
            moments, est = zarrnii_fibre_moments(znimg, psf_fwhm_um=None, **kwargs)
        print("moments computed")
        gammas = [float(g) for g in est["gammas"]]
        print(f"gammas {np.round(gammas, 4)}, threshold {est['threshold']:.4f}")
        if "psf_fwhm_um" in est:
            print(
                f"PSF FWHM (z, y, x) {np.round(est['psf_fwhm_um'], 2)} um, "
                f"equalised by extra sigma {np.round(est['psf_sigma_um'], 2)} um"
            )

        fig = make_qc_figure(
            est,
            scales_um,
            threshold_label(
                params.threshold_method,
                params.threshold_factor,
                params.threshold_quantile,
            ),
        )
        fig.savefig(snakemake.output.qc_png)

        with open(snakemake.output.params_json, "w") as f:
            json.dump(
                {
                    "scales_um": [float(s) for s in scales_um],
                    # None when the scales were given in um (fibre_scales_um)
                    "scale_factors": (
                        None
                        if params.scales_um
                        else [float(s) / min(spacing_um) for s in scales_um]
                    ),
                    "level": filter_level,
                    "output_level": output_level,
                    "spacing_um": [float(s) for s in spacing_um],
                    "block_factor": list(factor),
                    "block_voxels": int(np.prod(factor)),
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
                    "threshold_factor": (
                        float(params.threshold_factor)
                        if params.threshold_method == "quantile"
                        else None
                    ),
                    "threshold_quantile": (
                        float(params.threshold_quantile)
                        if params.threshold_method == "quantile"
                        else None
                    ),
                    "edge_min_weight": float(params.edge_min_weight),
                    "psf_fwhm_um": (
                        [float(v) for v in est["psf_fwhm_um"]]
                        if "psf_fwhm_um" in est
                        else None
                    ),
                    "psf_sigma_um": (
                        [float(v) for v in est["psf_sigma_um"]]
                        if "psf_sigma_um" in est
                        else None
                    ),
                    "psf_n_puncta": est.get("psf_n_puncta"),
                    "soft_threshold": float(params.soft_lo),
                    "local_norm_um": est.get("local_norm_um"),
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
