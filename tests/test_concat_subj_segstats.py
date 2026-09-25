from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


def _find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").exists():
            return candidate
    raise RuntimeError("Could not locate repository root from test path")


REPO_ROOT = _find_repo_root(Path(__file__).parent)
SCRIPT_PATH = REPO_ROOT / "spimquant/workflow/scripts/concat_subj_segstats.py"

spec = spec_from_file_location("concat_subj_segstats", SCRIPT_PATH)
concat_mod = module_from_spec(spec)
spec.loader.exec_module(concat_mod)


def test_load_tsvs_with_metadata_merges_regionpropstats_like_tables(tmp_path):
    participants_df = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-02"],
            "treatment": ["vehicle", "drug"],
            "sex": ["F", "M"],
        }
    )

    sub01_dir = tmp_path / "sub-01" / "tabular"
    sub02_dir = tmp_path / "sub-02" / "tabular"
    sub01_dir.mkdir(parents=True)
    sub02_dir.mkdir(parents=True)

    sub01_path = (
        sub01_dir
        / "sub-01_seg-ctx_from-ABAv3_stain-Abeta_level-3_desc-threshold_regionpropstats.tsv"
    )
    sub02_path = (
        sub02_dir
        / "sub-02_seg-ctx_from-ABAv3_stain-Abeta_level-3_desc-threshold_regionpropstats.tsv"
    )

    pd.DataFrame(
        {
            "index": [1, 2],
            "name": ["RegionA", "RegionB"],
            "count": [3, 5],
            "volume_mean": [10.0, 20.0],
        }
    ).to_csv(sub01_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "index": [1, 2],
            "name": ["RegionA", "RegionB"],
            "count": [4, 6],
            "volume_mean": [11.0, 21.0],
        }
    ).to_csv(sub02_path, sep="\t", index=False)

    combined = concat_mod.load_tsvs_with_metadata(
        [str(sub01_path), str(sub02_path)], participants_df
    )

    assert list(combined["participant_id"]) == ["sub-01", "sub-01", "sub-02", "sub-02"]
    assert list(combined["name"]) == ["RegionA", "RegionB", "RegionA", "RegionB"]
    assert list(combined["treatment"]) == ["vehicle", "vehicle", "drug", "drug"]
    assert "sex" in combined.columns


def test_main_reads_tsv_files_input_key(tmp_path):
    participants_path = tmp_path / "participants.tsv"
    pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-02"],
            "treatment": ["vehicle", "drug"],
        }
    ).to_csv(participants_path, sep="\t", index=False)

    sub01_dir = tmp_path / "sub-01" / "tabular"
    sub02_dir = tmp_path / "sub-02" / "tabular"
    sub01_dir.mkdir(parents=True)
    sub02_dir.mkdir(parents=True)

    sub01_path = sub01_dir / "sub-01_desc-threshold_regionpropstats.tsv"
    sub02_path = sub02_dir / "sub-02_desc-threshold_regionpropstats.tsv"
    output_path = tmp_path / "group" / "merged.tsv"

    pd.DataFrame({"index": [1], "name": ["RegionA"], "count": [3]}).to_csv(
        sub01_path, sep="\t", index=False
    )
    pd.DataFrame({"index": [1], "name": ["RegionA"], "count": [4]}).to_csv(
        sub02_path, sep="\t", index=False
    )

    concat_mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(
            participants_tsv=str(participants_path),
            tsv_files=[str(sub01_path), str(sub02_path)],
        ),
        output=SimpleNamespace(merged_tsv=str(output_path)),
    )

    concat_mod.main()

    merged = pd.read_csv(output_path, sep="\t")
    assert list(merged["participant_id"]) == ["sub-01", "sub-02"]
    assert list(merged["count"]) == [3, 4]
    assert list(merged["treatment"]) == ["vehicle", "drug"]


def test_main_accepts_legacy_segstats_tsvs_input_key(tmp_path):
    participants_path = tmp_path / "participants.tsv"
    pd.DataFrame({"participant_id": ["sub-01"], "treatment": ["vehicle"]}).to_csv(
        participants_path, sep="\t", index=False
    )

    sub01_dir = tmp_path / "sub-01" / "tabular"
    sub01_dir.mkdir(parents=True)
    sub01_path = sub01_dir / "sub-01_desc-threshold_mergedsegstats.tsv"
    output_path = tmp_path / "group" / "merged.tsv"

    pd.DataFrame({"index": [1], "name": ["RegionA"], "count": [3]}).to_csv(
        sub01_path, sep="\t", index=False
    )

    concat_mod.snakemake = SimpleNamespace(
        input=SimpleNamespace(
            participants_tsv=str(participants_path),
            segstats_tsvs=[str(sub01_path)],
        ),
        output=SimpleNamespace(merged_tsv=str(output_path)),
    )

    concat_mod.main()

    merged = pd.read_csv(output_path, sep="\t")
    assert list(merged["participant_id"]) == ["sub-01"]
    assert list(merged["count"]) == [3]
