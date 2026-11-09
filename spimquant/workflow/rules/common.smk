"""
Common utility functions for SPIMquant workflows.

This module provides shared helper functions used across multiple workflow rules:

- resources_path(): Path resolver for resources directory
- get_template_path(): Template file locator with optional cropping
- get_template_for_reg(): Registration-specific template selector
- get_stains_all_subjects(): Stain consistency validator across subjects
- spim_size_gb(), spim_chunk_mb(): Image and chunk size of a job's SPIM input, for resources
- scaled_mem_mb(), scaled_runtime(): Resources that grow with a job's SPIM input size
- bids_oz_out(): bids() wrapper for OME-Zarr outputs (adds directory() when needed)
- bids_oz_in(): bids() wrapper for OME-Zarr inputs

These functions abstract away path complexity and ensure consistent file
organization following BIDS conventions.
"""

import functools
import json
import math
from pathlib import Path

import numpy as np


def resources_path(path):
    """Get path relative to the resources folder.

    Remote URLs pass through untouched unless the file has been
    pre-downloaded into resources/imported/ (matching basename), in which
    case the local copy is used. This skips the storage plugin's remote
    existence check at DAG-build time, which fails on hosts without
    internet access or when zenodo rate-limits the burst of checks.
    """
    resources = Path(workflow.basedir).parent / "resources"
    if path.startswith(("http://", "https://")):
        local = resources / "imported" / Path(path).name
        if local.exists():
            return str(local)
        return path
    else:
        return str(resources / path)


def bids_oz_out(**kwargs):
    """Generate a bids path for an OME-Zarr output file.

    Substitutes the ``{ext}`` placeholder in the suffix with the configured
    ``zarr_store_ext`` value (``ome.zarr`` or ``ozx``).  When the extension is
    ``ome.zarr`` the path is wrapped with ``directory()`` as required by
    Snakemake for directory outputs; for single-file formats (e.g. ``ozx``)
    the plain path string is returned.
    """
    suffix = kwargs.pop("suffix", "")
    path = bids(**kwargs, suffix=suffix.replace("{ext}", zarr_store_ext))
    if zarr_store_ext == "ome.zarr":
        return directory(path)
    return path


def bids_oz_in(**kwargs):
    """Generate a bids path for an OME-Zarr input file.

    Substitutes the ``{ext}`` placeholder in the suffix with the configured
    ``zarr_store_ext`` value (``ome.zarr`` or ``ozx``).  No ``directory()``
    wrapping is applied because Snakemake does not require it for inputs.
    """
    suffix = kwargs.pop("suffix", "")
    return bids(**kwargs, suffix=suffix.replace("{ext}", zarr_store_ext))


def seg_mask_input(**kwargs):
    """Input function for the segmentation mask of ``{desc}``.

    An object-level desc (e.g. ``dog+s3r4k6+inst``, see ``--object_method``)
    maps to the mask it was derived from (``dog+s3r4k6``), so mask-based
    outputs exist under both descs without a second mask.
    """
    pattern = bids_oz_in(desc="{desc}", suffix="mask.{ext}", **kwargs)

    def _mask(wildcards):
        wc = dict(wildcards.items())
        wc["desc"] = mask_desc(wc["desc"])
        return pattern.format(**wc)

    return _mask


def get_template_path(root, template, template_crop=None, suffix=None):
    """Get template path, optionally cropped based on hemisphere"""
    if suffix is None:
        suffix = stain_for_reg if use_spim_template else "anat"

    if template_crop is not None:
        return bids(
            root=root,
            template=template,
            desc=f"{template_crop}crop",
            suffix=f"{suffix}.nii.gz",
        )
    else:
        return bids(root=root, template=template, suffix=f"{suffix}.nii.gz")


def get_template_for_reg(wildcards):
    """Get the appropriate template file for registration, cropped if specified"""
    return get_template_path(root, wildcards.template, config.get("template_crop"))


def get_stains_all_subjects(ignore_stains=None):
    """Get stains across all subjects, optionally filtering out ignored stains."""
    ignore_set = set(ignore_stains) if ignore_stains else set()

    stain_sets = []
    for zarr in spim_raw.expand():
        channels = set(get_spim_channels(zarr))
        # Remove ignored stains
        channels = channels - ignore_set
        stain_sets.append(channels)

    if not stain_sets:
        raise ValueError(
            "No SPIM inputs found. Check your filters or participant_label."
        )
    reference = stain_sets[0]
    if all(s == reference for s in stain_sets):
        return sorted(reference)
    else:
        raise ValueError(
            f"Stains across subjects are not consistent: {stain_sets}\n\n Consider using --ignore-stains to ignore stains that are not consistent across subjects."
        )


def _open_spim_level(spim_path, level):
    """(dims, shape, chunks, dtype) of a stored pyramid level, or None if the
    pyramid has fewer levels; read from metadata only, without any image data."""
    if spim_path.endswith(".ims"):
        from imaris_ims_zarr import ImsProcessSafeStore

        if level >= ImsProcessSafeStore(spim_path).ResolutionLevels:
            return None
        store = ImsProcessSafeStore(spim_path, ResolutionLevelLock=level)
        return ("tczyx", store.shape, store.chunks, np.dtype(store.dtype))

    import zarr

    if spim_path.endswith((".ozx", ".zip")):
        group = zarr.open_group(zarr.storage.ZipStore(spim_path, mode="r"), mode="r")
    else:
        group = zarr.open_group(spim_path, mode="r")
    attrs = group.attrs.asdict()
    multiscales = attrs.get("ome", attrs)["multiscales"][0]
    datasets = multiscales["datasets"]
    if level >= len(datasets):
        return None
    arr = group[datasets[level]["path"]]
    dims = [axis["name"] for axis in multiscales["axes"]]
    return (dims, arr.shape, arr.chunks, np.dtype(arr.dtype))


@functools.cache
def get_spim_info(spim_path, level=0):
    """Shape and chunk shape (by axis name) and bytes per voxel of a SPIM input.

    Read from the image metadata of .ome.zarr, .ozx, .ims or remote stores, where
    the file size on disk says little about the image size (compressed, chunked,
    or not local at all). Only metadata is read, so this is cheap enough for
    resource functions evaluated by the scheduler. Levels beyond the stored
    pyramid are halved in each spatial dimension per level, as ZarrNii does.
    """
    spim_path, level = str(spim_path), int(level)
    stored = level
    while (info := _open_spim_level(spim_path, stored)) is None and stored > 0:
        stored -= 1
    if info is None:
        raise ValueError(f"No pyramid levels found in {spim_path}")
    dims, shape, chunks, dtype = info
    factor = 2 ** (level - stored)
    spatial = {"x", "y", "z"}
    sizes = {
        d: ((max(1, n // factor), max(1, c // factor)) if d in spatial else (n, c))
        for d, n, c in zip(dims, shape, chunks)
        if d != "t"
    }
    # As in dask, a chunk is no larger than the image
    return {
        "shape": {d: n for d, (n, c) in sizes.items()},
        "chunks": {d: min(n, c) for d, (n, c) in sizes.items()},
        "itemsize": dtype.itemsize,
    }


def _spatial_bytes(sizes, itemsize):
    return math.prod(n for dim, n in sizes.items() if dim != "c") * itemsize


def get_spim_size_gb(spim_path, level=0):
    """Uncompressed size of one channel of a SPIM input at `level`, in GB.

    For block-parallel (dask) rules this mostly sets the run time."""
    info = get_spim_info(spim_path, level)
    return _spatial_bytes(info["shape"], info["itemsize"]) / 1e9


def get_spim_chunk_mb(spim_path, level=0):
    """Uncompressed size of one chunk of one channel of a SPIM input, in MB.

    For block-parallel (dask) rules this, times the threads, mostly sets the
    peak memory."""
    info = get_spim_info(spim_path, level)
    return _spatial_bytes(info["chunks"], info["itemsize"]) / 1e6


def spim_size_gb(wildcards, level=0):
    """get_spim_size_gb() for the SPIM input of a job, e.g. to scale its resources."""
    return get_spim_size_gb(spim_input(wildcards), level)


def spim_chunk_mb(wildcards, level=0):
    """get_spim_chunk_mb() for the SPIM input of a job, e.g. to scale its resources."""
    return get_spim_chunk_mb(spim_input(wildcards), level)


def _job_spim_size_gb(wildcards, level):
    """spim_size_gb() at `level`, or at the job's {level} wildcard if None."""
    return spim_size_gb(wildcards, wildcards.level if level is None else level)


def _cap_mem_mb(mem_mb):
    """mem_mb, at most config max_mem_mb (the memory of a node)."""
    cap = config.get("max_mem_mb")
    return int(min(mem_mb, cap) if cap else mem_mb)


def scaled_mem_mb(base_mb, mb_per_gb, level=None):
    """mem_mb for a rule whose peak memory grows with the size of its SPIM input.

    base_mb + mb_per_gb * GB of one channel of the job's SPIM input at `level`
    (default: the job's {level} wildcard), times the attempt, so that a job
    retried with --retries gets more memory, up to config max_mem_mb. Fitted
    from --benchmark runs, see tools/summarize_benchmarks.py."""

    def mem_mb(wildcards, attempt):
        return _cap_mem_mb(
            (base_mb + mb_per_gb * _job_spim_size_gb(wildcards, level)) * attempt
        )

    return mem_mb


def scaled_runtime(base_min, min_per_gb, level=None):
    """runtime (minutes) for a rule whose run time grows with the size of its
    SPIM input; as scaled_mem_mb()."""

    def runtime(wildcards, attempt):
        return int(
            (base_min + min_per_gb * _job_spim_size_gb(wildcards, level)) * attempt
        )

    return runtime


def reg_mem_mb(base):
    """mem_mb for a registration rule: base, doubled for a multichannel method, times
    the attempt, so that a job retried with --retries gets more."""

    def mem_mb(wildcards, attempt):
        multichannel = reg_presets[wildcards.reg].get("multichannel")
        return _cap_mem_mb(base * (2 if multichannel else 1) * attempt)

    return mem_mb


def per_attempt(base):
    """A resource for a rule whose needs don't follow its SPIM input size: base
    times the attempt, so that a job retried with --retries gets more."""

    def resource(wildcards, attempt):
        return base * attempt

    return resource


def get_spim_json_path(spim_path):
    spim_path = str(spim_path)
    if spim_path.endswith(".ozx"):
        return spim_path[: -len(".ozx")] + ".json"
    if spim_path.endswith(".ome.zarr.zip"):
        return spim_path[: -len(".ome.zarr.zip")] + ".json"
    if spim_path.endswith(".ome.zarr"):
        return spim_path[: -len(".ome.zarr")] + ".json"
    return str(Path(spim_path).with_suffix(".json"))


def spim_json_input(wildcards):
    """JSON sidecar of a job's SPIM input, or nothing if it has none.

    Sidecar values (e.g. orientation) change the outputs of rules that read the
    SPIM input, so any rule with spim_input as an input should have this as well,
    to rerun if the sidecar changes."""
    json_path = get_spim_json_path(spim_input(wildcards))
    return json_path if Path(json_path).exists() else []


def _normalize_channel_labels(value):
    if value is None:
        return None
    if isinstance(value, str):
        return [value]
    if isinstance(value, list):
        channels = [str(item) for item in value if isinstance(item, (str, int, float))]
        return channels or None
    return None


def _extract_channels_from_json(metadata):
    """Extract channel labels from sidecar keys, using "SampleStaining" key"""
    channels = _normalize_channel_labels(metadata.get("SampleStaining"))
    if channels:
        return channels

    return None


def get_spim_json_overrides(spim_path):
    json_path = get_spim_json_path(spim_path)
    json_file = Path(json_path)
    if not json_file.exists():
        return {}

    with json_file.open() as fp:
        metadata = json.load(fp)

    overrides = {}
    channels = _extract_channels_from_json(metadata)
    if channels:
        overrides["set_channel_labels"] = channels

    orientation = metadata.get("OrientationStringXYZ")
    if orientation:
        overrides["orientation"] = str(orientation)

    return overrides


def get_spim_channels(spim_path):
    overrides = get_spim_json_overrides(spim_path)
    if "set_channel_labels" in overrides:
        return overrides["set_channel_labels"]

    return ZarrNii.from_file(spim_path).list_channels()


def get_mri_reg_targets(pattern, **kwargs):
    """Expand a target pattern for SPIM subjects that also have MRI data.

    Uses ``inputs['mri']`` to determine which subjects have MRI data, then
    filters the reference scans to those subjects before expanding *pattern*.
    This prevents ``all_mri_reg`` from generating targets for SPIM subjects
    that lack a corresponding MRI, which would cause downstream rules (e.g.
    ``average_mri``) to fail with empty input lists.

    Parameters
    ----------
    pattern:
        A bids path string containing wildcards.
    **kwargs:
        Extra keyword arguments forwarded to :meth:`BidsComponent.expand`
        (e.g. ``template=config["template"]``).
    """
    if "mri" not in inputs:
        return []
    mri_subjects = set(inputs["mri"].entities.get("subject", []))
    # additional scans (--aux_spim) are not registered to the MRI, only reference scans
    return ref_spim.filter(subject=mri_subjects).expand(pattern, **kwargs)


def get_regionprops_parquet(wildcards):
    if do_vessels:
        return bids(
            root=root,
            datatype="tabular",
            desc="{desc}+vessels",
            space="{template}",
            vessels=stain_for_vessels,
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        )
    else:
        return bids(
            root=root,
            datatype="tabular",
            desc="{desc}",
            space="{template}",
            suffix="regionprops.parquet",
            **inputs["spim"].wildcards,
        )


# --- registration methods: named presets, used as the desc- of their outputs ---


def _resolve_reg_preset(name):
    presets = config["registration"]["presets"]
    if name not in presets:
        raise ValueError(
            f"Unknown registration method {repr(name)}, choose from {sorted(presets)}"
        )
    preset = dict(presets[name])
    base = preset.pop("base", None)
    return {**_resolve_reg_preset(base), **preset} if base else preset


reg_presets = {
    name: _resolve_reg_preset(name) for name in config["registration"]["presets"]
}
for _name, _preset in reg_presets.items():
    if not _name.isalnum() or _name == "cohort":
        raise ValueError(
            f"Registration method name {repr(_name)} must be alphanumeric, and "
            "cohort is reserved for warps from a cohort template build"
        )
    if _preset["backend"] not in ("greedy", "ants", "fireants"):
        raise ValueError(f"Unknown backend in registration method {_name}: {_preset}")
# every method has a multichannel variant, named <method>Multi
for _name, _preset in list(reg_presets.items()):
    if not _preset.get("multichannel"):
        reg_presets.setdefault(f"{_name}Multi", {**_preset, "multichannel": True})


def select_registration_methods(reg_stains):
    """Set the subject registration methods, given the stains shared with the template.

    With two or more, --registration_method becomes its multichannel variant, unless
    --registration_single_channel. With fewer, multichannel methods can't be used.
    """
    method = config["registration_method"]
    if method not in reg_presets:
        raise ValueError(
            f"Unknown registration method {repr(method)}, choose from "
            f"{sorted(reg_presets)}"
        )
    if len(reg_stains) > 1:
        if config["registration_single_channel"]:
            print(f"Registering to the template with one stain: {reg_stains[0]}")
        else:
            if not reg_presets[method].get("multichannel"):
                config["registration_method"] = f"{method}Multi"
            print(
                f"Registering with {config['registration_method']}, using the "
                f"stains shared with the template: {reg_stains}"
            )
        return
    if reg_presets[method].get("multichannel"):
        raise ValueError(
            f"Registration method {method} is multichannel, but the template shares "
            f"fewer than two stains with the subjects ({reg_stains}). Choose a "
            "single-channel method, or a template with SPIM templates of their stains."
        )
    config["registration_tune_methods"] = [
        m
        for m in config["registration_tune_methods"]
        if not reg_presets[m].get("multichannel")
    ]


def reg_constraint(backend=None):
    """Regex matching the registration method names, optionally of one backend."""
    names = [n for n, p in reg_presets.items() if backend in (None, p["backend"])]
    return "|".join(names) if names else "(?!)"


def reg_iters(preset, key):
    """Iterations of a preset, or very few when running with --sloppy."""
    if not config["sloppy"]:
        return preset[key]
    return [10, 0, 0] if preset["backend"] == "fireants" else "10x0x0"


def reg_note(preset):
    """One-line summary of a preset for reports."""
    return ", ".join(f"{k}: {v}" for k, v in preset.items())


def ants_reg_args(preset, fixed, moving, center_init, extra_fixed=(), extra_moving=()):
    """antsRegistration stages for a preset (the output prefix is added by the rule).

    Each image in extra_fixed adds a metric term, against the matching image in
    extra_moving, or against the same moving image if there is no extra_moving.
    """
    pairs = [
        (fixed, moving),
        *zip(extra_fixed, extra_moving or [moving] * len(extra_fixed)),
    ]

    def metric_args(metric):
        forms = {
            "MI": "MI[{},{},1,32,Regular,0.25]",
            "CC": "CC[{},{},1," + str(preset["cc_radius"]) + "]",
        }
        return " ".join("-m " + forms[metric].format(*pair) for pair in pairs)

    levels = f"-f {preset['shrink']} -s {preset['smooth']}vox "
    args = "-d 3 --float 1 --collapse-output-transforms 1 -v 1 -n Linear "
    if center_init:
        args += f"--initial-moving-transform [{fixed},{moving},0] "  # image centers
    if preset["affine"]:
        args += (
            f"-t Affine[0.1] {metric_args(preset['affine_metric'])} "
            f"-c [{reg_iters(preset, 'affine_iters')},1e-6,10] {levels}"
        )
    return args + (
        f"-t SyN[{preset['syn']}] {metric_args(preset['metric'])} "
        f"-c [{reg_iters(preset, 'deform_iters')},1e-6,10] {levels}"
    )


# printf commands for the identity transforms of presets without an affine stage
reg_identity_ras = "printf '1 0 0 0\\n0 1 0 0\\n0 0 1 0\\n0 0 0 1\\n'"
reg_identity_itk = (
    "printf '#Insight Transform File V1.0\\n#Transform 0\\n"
    "Transform: MatrixOffsetTransformBase_double_3_3\\n"
    "Parameters: 1 0 0 0 1 0 0 0 1 0 0 0\\nFixedParameters: 0 0 0\\n'"
)


# --- subject warps from a cohort template build, when --template is that cohort ---

cohort_warps = config["templates"][config["template"]].get("cohort_warps")


def _cohort_warp_paths(wildcards):
    """Files that make up the subject<->cohort template warps of a subject."""
    entities = {k: wildcards[k] for k in inputs["spim"].wildcards}
    row = os.path.basename(bids(stain="x", suffix="x", **entities)).split("_stain-")[0]
    stem = os.path.join(
        cohort_warps["warps_dir"], f"{row}_stain-{cohort_warps['reg_stain']}"
    )
    orig = cohort_warps["orig_template"]
    return {
        "warp": f"{stem}_warp.nii.gz",
        "invwarp": f"{stem}_invwarp.nii.gz",
        "to_orig": bids(
            root=cohort_warps["orig_root"],
            datatype="xfm",
            from_="subject",
            to=orig,
            desc="affine",
            suffix="xfm.nii.gz",
            **entities,
        ),
        "from_orig": bids(
            root=cohort_warps["orig_root"],
            datatype="xfm",
            from_=orig,
            to="subject",
            desc="affine",
            suffix="xfm.nii.gz",
            **entities,
        ),
    }


def _scan_key(wildcards):
    return tuple(getattr(wildcards, w) for w in _spim_key_wildcards)


def is_aux_scan(wildcards):
    """Whether a job's scan is an additional scan (--aux_spim)."""
    return _scan_key(wildcards) in aux_ref_keys


def ref_scan_wildcards(wildcards):
    """Wildcards of a job with its additional scan replaced by its reference scan."""
    ref = dict(zip(_spim_key_wildcards, aux_ref_keys[_scan_key(wildcards)]))
    return {**dict(wildcards.items()), **ref}


def ref_scan_input(pattern):
    """Input function for `pattern` of the reference scan of an additional scan."""

    def _input(wildcards):
        return pattern.format(**ref_scan_wildcards(wildcards))

    return _input


def subject_reg_method(wildcards):
    """Registration desc for a subject: its cohort build warps if any, else the method.

    Additional scans (--aux_spim) get theirs through their reference scan (viaref).
    """
    if is_aux_scan(wildcards):
        return "viaref"
    if wildcards.template == config["template"] and cohort_warps:
        if all(os.path.exists(f) for f in _cohort_warp_paths(wildcards).values()):
            return "cohort"
    return config["registration_method"]
