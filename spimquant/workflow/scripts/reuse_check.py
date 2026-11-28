"""Settings checks for --reuse_derivatives (imported by rules/reuse.smk).

The BIDS path of a segmentation or vessel output names its stain, level and
method, but not every setting that went into it (e.g. --correction_method for
a gmm mask, --vessel_threshold for the vesselFM mask). Before reusing such
outputs, the settings below are compared between the current config and the
config the prior dataset was made with (its config/snakebids.yml, which holds
the resolved config of the latest run on that dataset). If any differ, or the
prior config lacks one, that group of outputs is not reused.

Keys with a dot are nested (``models.vesselfm`` is config["models"]["vesselfm"]).
"""

import os

import yaml

# intensity correction chain read by the threshold methods
_CORRECTION = [
    "correction_method",
    "correction_level",
    "correction_n4_iters",
    "correction_n4_spline_spacing",
]

_LANTERN_INFERENCE = [
    "lantern_config",
    "models.lantern",
    "plaque_level",
    "plaque_iso_res",
    "plaque_n_folds",
    "plaque_tile",
    "plaque_stride",
    "plaque_merge_mode",
]

# group -> settings that change its outputs but not their paths
REUSE_SETTINGS = {
    # th<N>, otsu+k<k>i<i>, gmm+n<n>k<k> masks
    "threshold": _CORRECTION + ["seg_hist_bin_width", "seg_hist_percentile_range"],
    # dog+s<s>r<r>k<k> masks
    "dog": _CORRECTION
    + [
        "dog_pool_level_offset",
        "dog_pool_factor",
        "dog_pool_min_volume_um3",
        "dog_pool_lowres_dilate_um",
        "dog_pool_dilate_um",
        "dog_stat_blocks",
        "dog_stat_samples_per_block",
    ],
    # LANTERN vote map (probseg)
    "lantern_probseg": _LANTERN_INFERENCE,
    # LANTERN binary mask, made from the vote map
    "lantern_mask": _LANTERN_INFERENCE
    + [
        "plaque_vote_threshold",
        "plaque_fp_min_votes",
        "min_size",
        "label_filter",
        "label_filter_size",
    ],
    # vesselFM mask and everything made from it
    "vessels": [
        "models.vesselfm",
        "vessel_tile",
        "vessel_stride",
        "vessel_merge_mode",
        "vessel_threshold",
    ],
}

MISSING = object()


def get_setting(cfg, key):
    """Value of a (dotted) key in a config dict, or MISSING."""
    value = cfg
    for part in key.split("."):
        if not isinstance(value, dict) or part not in value:
            return MISSING
        value = value[part]
    return value


def _normalize(value):
    """Make YAML-loaded and in-memory values comparable (1e6 == 1000000.0)."""
    if isinstance(value, bool):
        return ("bool", value)  # True == 1.0 otherwise
    if value is None or value is MISSING:
        return value
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, (list, tuple)):
        return [_normalize(v) for v in value]
    if isinstance(value, dict):
        return {k: _normalize(v) for k, v in value.items()}
    return value


def setting_differences(prior_cfg, current_cfg, keys):
    """[(key, prior value, current value)] for keys that differ.

    A key missing from the prior config counts as a difference (MISSING), since
    the prior run's value cannot be known.
    """
    diffs = []
    for key in keys:
        prior = get_setting(prior_cfg, key)
        current = get_setting(current_cfg, key)
        if _normalize(prior) != _normalize(current):
            diffs.append((key, prior, current))
    return diffs


def format_differences(diffs):
    """Human-readable list of setting differences."""

    def show(value):
        return "(not recorded)" if value is MISSING else repr(value)

    return "; ".join(
        f"{key}: prior {show(prior)}, now {show(current)}"
        for key, prior, current in diffs
    )


def load_prior_config(reuse_root):
    """The resolved config of the prior dataset, or None if it has none."""
    path = os.path.join(reuse_root, "config", "snakebids.yml")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return yaml.safe_load(f) or {}
