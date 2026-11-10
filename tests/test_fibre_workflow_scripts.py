"""End-to-end test of the fibre_moments and fibre_metrics workflow scripts."""

import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

zarrnii = pytest.importorskip("zarrnii")
nib = pytest.importorskip("nibabel")


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


ROOT = _find_repo_root(Path(__file__).parent)
SCRIPTS_DIR = ROOT / "spimquant/workflow/scripts"


def _load(name):
    spec = spec_from_file_location(name, SCRIPTS_DIR / f"{name}.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def scripts():
    sys.path.insert(0, str(SCRIPTS_DIR))
    try:
        yield SimpleNamespace(
            moments=_load("fibre_moments"),
            metrics=_load("fibre_metrics"),
            roi=_load("fibre_roi_stats"),
            merge=_load("merge_fibre_roi_stats"),
            template=_load("fibre_metrics_template"),
            lib=_load("fibre_orientation"),
        )
    finally:
        sys.path.remove(str(SCRIPTS_DIR))


def _phantom(shape, spacing_um, seed=0):
    """Fibres along array axis 2 in one half, along axis 1 in the other,
    inside an ellipsoidal 'brain' (zero outside, as the corrected image)."""
    rng = np.random.default_rng(seed)
    grids = np.meshgrid(
        *[np.arange(n) * s for n, s in zip(shape, spacing_um)], indexing="ij"
    )
    coords = np.stack(grids, axis=-1)
    ext = np.asarray(shape) * np.asarray(spacing_um)
    img = np.zeros(shape, np.float32)
    for k in range(60):
        axis = 2 if k % 2 == 0 else 1
        d = np.zeros(3)
        d[axis] = 1.0
        p = rng.uniform(0.1, 0.9, 3) * ext
        r = coords - p
        perp = r - (r @ d)[..., None] * d
        line = 100 * np.exp(-np.sum(perp**2, -1) / (2 * 6.0**2))
        half = (
            coords[..., 0] < ext[0] / 2 if axis == 2 else coords[..., 0] >= ext[0] / 2
        )
        img = np.maximum(img, line * half)
    img += 20 + rng.normal(0, 2, shape)
    centre = ext / 2
    inside = np.sum(((coords - centre) / (0.48 * ext)) ** 2, -1) <= 1
    return np.where(inside, img, 0).astype(np.float32)


def test_moments_and_metrics(scripts, tmp_path):
    import dask.array as da

    moments_mod, metrics_mod = scripts.moments, scripts.metrics
    spacing_um = (4.0, 4.0, 4.0)
    shape = (64, 64, 64)
    img = _phantom(shape, spacing_um)
    corrected = tmp_path / "corrected.ome.zarr"
    zarrnii.ZarrNii.from_darr(
        da.from_array(img[None], chunks=(1, 32, 32, 32)),
        spacing=tuple(s / 1000 for s in spacing_um),
        origin=(1.0, 2.0, 3.0),
        axes_units={"x": "millimeter", "y": "millimeter", "z": "millimeter"},
    ).to_ome_zarr(str(corrected), max_layer=3)

    out = SimpleNamespace(
        moments=str(tmp_path / "moments.ome.zarr"),
        params_json=str(tmp_path / "params.json"),
        qc_png=str(tmp_path / "qc.png"),
    )
    moments_mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(corrected=str(corrected)),
        output=out,
        params=SimpleNamespace(
            scales_um=[6, 10],
            block_um=32,
            lmax=4,
            alpha=0.001,
            beta=1.0,
            gamma_quantile=0.999,
            gamma_factor=0.5,
            threshold_method="yen",
            stat_blocks=4,
            stat_samples_per_block=5000,
            stat_level_offset=2,
            edge_min_weight=0.5,
            chunk_size=32,
            slab_chunks=1,
            dask_threads_per_worker=1,
            dask_reserve_mb=0,
        ),
        threads=2,
        resources={},
        config={"dask_scheduler": "threads", "zarrnii_out_kwargs": {}},
    )
    moments_mod.main()

    params = json.loads(Path(out.params_json).read_text())
    assert params["block_factor"] == [8, 8, 8]
    assert len(params["gammas"]) == 2
    assert Path(out.qc_png).exists()

    names = list(metrics_mod.SCALAR_METRICS) + ["dir", "dec", "odf", "moments_nii"]
    mout = SimpleNamespace(**{n: str(tmp_path / f"{n}.nii.gz") for n in names})
    metrics_mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(moments=out.moments), output=mout, threads=2
    )
    metrics_mod.main()

    ref = nib.load(mout.fa)
    assert ref.shape == (8, 8, 8)
    # block grid: 8x the voxel size (mm), centred on the blocks
    np.testing.assert_allclose(np.abs(np.diag(ref.affine)[:3]), 0.032, rtol=1e-5)
    dens = np.asarray(nib.load(mout.fibredensity).dataobj)
    fa = np.asarray(ref.dataobj)
    assert dens.max() > 0.01
    # NIfTI is x, y, z: array axis 0 (zarr z) is NIfTI axis 2. Fibres along
    # array axis 2 (world x) in the lower half of axis 0, along axis 1
    # (world y) in the upper half.
    direction = np.asarray(nib.load(mout.dir).dataobj)
    assert direction.shape == (8, 8, 8, 3)
    has = (dens > 0.02) & (fa > 0.5)
    lower = has.copy()
    lower[:, :, 4:] = False
    upper = has.copy()
    upper[:, :, :4] = False
    assert lower.sum() > 5 and upper.sum() > 5
    assert np.median(np.abs(direction[..., 0][lower])) > 0.9
    assert np.median(np.abs(direction[..., 1][upper])) > 0.9
    odf = nib.load(mout.odf)
    assert odf.shape == (8, 8, 8, 15)

    # ---- atlas-region stats: two labels splitting the volume along axis 0
    import pandas as pd

    moments = np.asarray(
        zarrnii.ZarrNii.from_file(out.moments).data.compute(), dtype=np.float64
    )
    labels = np.ones(moments.shape[1:], dtype=np.float32)
    labels[4:] = 2
    template_zn = zarrnii.ZarrNii.from_file(out.moments)
    dseg_path = str(tmp_path / "dseg.nii.gz")
    metrics_mod.write_nifti(template_zn, labels, dseg_path)
    label_tsv = tmp_path / "dseg.tsv"
    pd.DataFrame(
        {"index": [0, 1, 2, 3], "name": ["bg", "lower", "upper", "absent"]}
    ).to_csv(label_tsv, sep="\t", index=False)
    roi_tsv = tmp_path / "fibrestats.tsv"
    scripts.roi.snakemake = SimpleNamespace(
        input=SimpleNamespace(
            moments=out.moments, dseg=dseg_path, label_tsv=str(label_tsv)
        ),
        output=SimpleNamespace(tsv=str(roi_tsv)),
        params=SimpleNamespace(metrics=list(metrics_mod.SCALAR_METRICS)),
        threads=1,
    )
    scripts.roi.main()
    roi = pd.read_csv(roi_tsv, sep="\t").set_index("name")
    # pooled moments: identical to summing the moments by hand
    for lab, name in ((1, "lower"), (2, "upper")):
        pooled = moments[:, labels.astype(int) == lab].sum(axis=1)
        expected = scripts.lib.moments_to_metrics(pooled)
        for m in metrics_mod.SCALAR_METRICS:
            assert roi.loc[name, m] == pytest.approx(float(expected[m]), abs=1e-5)
    # one fibre direction per half: well aligned within each region
    assert roi.loc["lower", "odi"] < 0.1 and roi.loc["upper", "odi"] < 0.1
    assert roi.loc["absent", "fibrevoxels"] == 0
    assert roi.loc["absent", "odi"] == 0

    # ---- merged across stains, segstats-style column names
    merged = scripts.merge.merge_tables(
        [pd.read_csv(roi_tsv, sep="\t")] * 2, ["VAChT", "MOG"]
    )
    assert "VAChT+odi" in merged.columns and "MOG+fibredensity" in merged.columns
    assert len(merged) == 4

    # ---- template-space metrics from the moments NIfTI (identity warp)
    tpl_metrics = [m for m in metrics_mod.SCALAR_METRICS if m != "gfa"]
    tout = SimpleNamespace(
        **{m: str(tmp_path / f"tpl_{m}.nii.gz") for m in tpl_metrics}
    )
    scripts.template.snakemake = SimpleNamespace(
        input=SimpleNamespace(moments_nii=mout.moments_nii),
        output=tout,
        params=SimpleNamespace(metrics=tpl_metrics, min_tissue_fraction=0.0),
    )
    scripts.template.main()
    assert nib.load(mout.moments_nii).shape[-1] == scripts.lib.SH_START
    for m in ("fibredensity", "frangi", "scale", "odi", "odipri", "fa"):
        a = nib.load(getattr(tout, m))
        b = nib.load(getattr(mout, m))
        np.testing.assert_allclose(a.affine, b.affine, atol=1e-9)
        np.testing.assert_allclose(
            np.asarray(a.dataobj), np.squeeze(np.asarray(b.dataobj)), atol=1e-5
        )


def test_slabwise_equals_whole(scripts):
    """Slab-wise reading gives the same maps and pooled moments."""
    da = pytest.importorskip("dask.array")
    lib = scripts.lib
    rng = np.random.default_rng(1)
    m = np.abs(rng.normal(size=(lib.n_moments(2), 7, 5, 6))) + 1
    darr = da.from_array(m, chunks=(4, 3, 5, 6))
    slabs = list(scripts.metrics.z_slabs(darr, slab_bytes=8 * m.shape[0] * 5 * 6 * 2))
    assert len(slabs) == 4  # 2 z-planes per slab
    whole = scripts.metrics.compute_maps(m)
    for z0, z1, part in slabs:
        res = scripts.metrics.compute_maps(part)
        for k in ("odi", "fa", "gfa", "odf"):
            np.testing.assert_allclose(res[k], whole[k][..., z0:z1, :, :], rtol=1e-5)
    labels = rng.integers(0, 4, size=m.shape[1:])
    pooled = scripts.roi.pooled_moments(darr, labels, 4, slab_bytes=8 * m[:, :1].size)
    np.testing.assert_allclose(pooled, lib.aggregate_moments(m, labels, 4))


def test_block_factor_snaps_to_chunk(scripts):
    bf = scripts.moments.block_factor
    # level 1 of a 4x mouse brain: 30 um -> 7.5 and 9.2 voxels
    assert bf(30, (4.0, 3.25, 3.25)) == (8, 9, 9)
    assert bf(30, (4.0, 3.25, 3.25), 128) == (8, 8, 8)
    assert bf(100, (4.0, 1.0, 1.0), 128) == (32, 128, 128)
