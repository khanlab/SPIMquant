"""Settings checks for reusing segmentation/vessel outputs (reuse_check.py)."""

import sys
from pathlib import Path

import yaml


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


REPO = _find_repo_root(Path(__file__).parent)
sys.path.insert(0, str(REPO / "spimquant/workflow/scripts"))

from reuse_check import (  # noqa: E402
    MISSING,
    REUSE_SETTINGS,
    format_differences,
    get_setting,
    load_prior_config,
    setting_differences,
)


def test_get_setting_nested_and_missing():
    cfg = {"models": {"vesselfm": "url"}, "vessel_threshold": 0.5}
    assert get_setting(cfg, "models.vesselfm") == "url"
    assert get_setting(cfg, "vessel_threshold") == 0.5
    assert get_setting(cfg, "models.lantern") is MISSING
    assert get_setting(cfg, "vessel_threshold.x") is MISSING


def test_equal_settings_after_yaml_round_trip():
    current = {
        "dog_pool_min_volume_um3": 1000000.0,
        "seg_hist_percentile_range": [0.5, 99.9],
        "correction_level": 6,
        "label_filter": False,
        "models": {"vesselfm": "url"},
    }
    prior = yaml.safe_load(
        "dog_pool_min_volume_um3: 1.0e+6\n"
        "seg_hist_percentile_range: [0.5, 99.9]\n"
        "correction_level: 6.0\n"
        "label_filter: false\n"
        "models: {vesselfm: url}\n"
    )
    assert (
        setting_differences(prior, current, list(current) + ["models.vesselfm"]) == []
    )


def test_differences_and_missing_prior_key():
    current = {"correction_method": "n4", "correction_level": 6, "min_size": 0}
    prior = {"correction_method": "gaussian", "correction_level": 6}
    diffs = setting_differences(prior, current, list(current))
    assert diffs == [
        ("correction_method", "gaussian", "n4"),
        ("min_size", MISSING, 0),
    ]
    text = format_differences(diffs)
    assert "correction_method: prior 'gaussian', now 'n4'" in text
    assert "min_size: prior (not recorded), now 0" in text


def test_bool_is_not_number():
    assert setting_differences(
        {"label_filter": 1}, {"label_filter": True}, ["label_filter"]
    )


def test_mask_settings_include_vote_map_settings():
    # a changed vote map setting must also stop reuse of the mask made from it
    assert set(REUSE_SETTINGS["lantern_probseg"]) < set(REUSE_SETTINGS["lantern_mask"])


def test_settings_are_config_keys():
    """Every checked setting exists in the default config."""
    cfg = yaml.safe_load((REPO / "spimquant/config/snakebids.yml").read_text())
    for key in {k for keys in REUSE_SETTINGS.values() for k in keys}:
        top = key.split(".")[0]
        assert get_setting(cfg, key) is not MISSING or any(
            opt.get("dest", name.lstrip("-").replace("-", "_")) == top
            for name, opt in cfg["parse_args"].items()
        ), key


def test_load_prior_config(tmp_path):
    assert load_prior_config(str(tmp_path)) is None
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "snakebids.yml").write_text("vessel_threshold: 0.5\n")
    assert load_prior_config(str(tmp_path)) == {"vessel_threshold": 0.5}
