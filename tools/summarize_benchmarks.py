#!/usr/bin/env python
"""Summarize SPIMquant job benchmarks, for tuning the resources of each rule.

Run SPIMquant with --benchmark, then point this at one or more output folders
(e.g. one per dataset size). Every job's benchmark is collected into a per-job
table, and a per-rule table compares what was used against what was requested:

  mem_mb_req / max_rss_mb      requested vs peak memory (max over jobs)
  runtime_req / max_min        requested vs longest wall time, in minutes
  threads / cpu_eff            requested threads, and how many were kept busy
                               (cpu_time / (wall time * threads))

  spim_gb_l0 / spim_gb         size of one channel of the job's SPIM input at
                               level 0, and at the job's {level} wildcard if it
                               has one (from benchmarks/spim_sizes.tsv)
  chunk_mb                     size of one chunk of one channel of that input,
                               at the job's {level} (or level 0)

For block-parallel (dask) rules, peak memory should follow threads * chunk_mb
and run time spim_gb. Rules reading SPIMquant's own zarr outputs see the chunks
in benchmarks/run_config.json (zarr_store_chunks) instead of the input's.

The peak memory of a job includes all of its child processes (e.g. dask workers).
Run times from a shared node are noisy; memory numbers are not.

Usage:
  python tools/summarize_benchmarks.py OUT_DIR [OUT_DIR ...] [-o PREFIX]

Writes PREFIX_jobs.tsv and PREFIX_rules.tsv (default PREFIX: benchmarks) and
prints the per-rule table.
"""

import argparse
import json
from pathlib import Path

import pandas as pd


def load_spim_sizes(out_dir):
    """{(scan wildcards, level): (GB, chunk MB) per channel}, scan wildcard names."""
    path = Path(out_dir, "benchmarks", "spim_sizes.tsv")
    if not path.exists():
        return {}, []
    # As strings, so that e.g. subject 01 is not read as 1
    sizes = pd.read_csv(path, sep="\t", dtype=str)
    key_cols = list(sizes.columns[: sizes.columns.get_loc("level")])
    lookup = {
        (tuple(row[key_cols]), int(row["level"])): (
            float(row["gb_per_channel"]),
            float(row["chunk_mb_per_channel"]),
        )
        for _, row in sizes.iterrows()
    }
    return lookup, key_cols


def get_spim_size(wildcards, sizes, key_cols, level):
    """(GB, chunk MB) of the job's SPIM input at `level`, or (None, None)."""
    if not key_cols or not all(c in wildcards for c in key_cols):
        return None, None
    key = tuple(str(wildcards[c]) for c in key_cols)
    return sizes.get((key, int(level)), (None, None))


def load_jobs(out_dir, label):
    sizes, key_cols = load_spim_sizes(out_dir)
    rows = []
    for path in sorted(Path(out_dir, "benchmarks").rglob("*.jsonl")):
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            resources = rec.get("resources") or {}
            input_sizes = rec.get("input_size_mb") or {}
            wildcards = rec.get("wildcards") or {}
            level = str(wildcards.get("level", ""))
            spim_gb_l0, chunk_mb_l0 = get_spim_size(wildcards, sizes, key_cols, 0)
            spim_gb, chunk_mb = (
                get_spim_size(wildcards, sizes, key_cols, level)
                if level.isdigit()
                else (None, chunk_mb_l0)
            )
            rows.append(
                {
                    "dataset": label,
                    "rule": rec.get("rule_name"),
                    "wildcards": json.dumps(wildcards),
                    "s": rec.get("s"),
                    "max_rss_mb": rec.get("max_rss"),
                    "max_pss_mb": rec.get("max_pss"),
                    "max_uss_mb": rec.get("max_uss"),
                    "cpu_time_s": rec.get("cpu_time"),
                    "io_in_mb": rec.get("io_in"),
                    "io_out_mb": rec.get("io_out"),
                    "threads": rec.get("threads"),
                    "mem_mb_req": resources.get("mem_mb"),
                    "runtime_req": resources.get("runtime"),
                    "gpus_req": resources.get("gpus", resources.get("gpu")),
                    # Directories (e.g. .ome.zarr stores) are not summed, only stat'd
                    "input_mb": sum(input_sizes.values()),
                    "spim_gb_l0": spim_gb_l0,
                    "spim_gb": spim_gb,
                    "chunk_mb": chunk_mb,
                    "benchmark": str(path),
                }
            )
    return rows


def summarize_rules(jobs):
    jobs = jobs.copy()
    jobs["min"] = jobs["s"] / 60
    jobs["cpu_eff"] = jobs["cpu_time_s"] / (jobs["s"] * jobs["threads"])
    rules = (
        jobs.groupby(["dataset", "rule"])
        .agg(
            n=("s", "size"),
            threads=("threads", "max"),
            cpu_eff=("cpu_eff", "median"),
            mem_mb_req=("mem_mb_req", "max"),
            max_rss_mb=("max_rss_mb", "max"),
            median_rss_mb=("max_rss_mb", "median"),
            runtime_req=("runtime_req", "max"),
            max_min=("min", "max"),
            median_min=("min", "median"),
            max_input_mb=("input_mb", "max"),
            spim_gb_l0=("spim_gb_l0", "max"),
            spim_gb=("spim_gb", "max"),
            chunk_mb=("chunk_mb", "max"),
        )
        .reset_index()
    )
    rules["mem_ratio"] = rules["mem_mb_req"] / rules["max_rss_mb"]
    rules["runtime_ratio"] = rules["runtime_req"] / rules["max_min"]
    return rules.sort_values(["dataset", "mem_mb_req"], ascending=[True, False])


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "out_dirs", nargs="+", help="SPIMquant output folders run with --benchmark"
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        help="Dataset label for each output folder (default: folder name)",
    )
    parser.add_argument(
        "-o", "--prefix", default="benchmarks", help="Prefix for the output tables"
    )
    args = parser.parse_args()

    labels = args.labels or [Path(d).resolve().name for d in args.out_dirs]
    if len(labels) != len(args.out_dirs):
        parser.error("--labels needs one label per output folder")

    rows = []
    for out_dir, label in zip(args.out_dirs, labels):
        rows += load_jobs(out_dir, label)
    if not rows:
        parser.error(
            "no benchmarks/**/*.jsonl found; was SPIMquant run with --benchmark?"
        )

    jobs = pd.DataFrame(rows)
    num_cols = [
        "s",
        "max_rss_mb",
        "max_pss_mb",
        "max_uss_mb",
        "cpu_time_s",
        "io_in_mb",
        "io_out_mb",
        "threads",
        "mem_mb_req",
        "runtime_req",
        "spim_gb_l0",
        "spim_gb",
        "chunk_mb",
    ]
    # Jobs too short for psutil to sample have "NA" for memory and CPU
    jobs[num_cols] = jobs[num_cols].apply(pd.to_numeric, errors="coerce")
    rules = summarize_rules(jobs)

    jobs.to_csv(f"{args.prefix}_jobs.tsv", sep="\t", index=False)
    rules.to_csv(f"{args.prefix}_rules.tsv", sep="\t", index=False)

    with pd.option_context(
        "display.max_rows",
        None,
        "display.width",
        250,
        "display.float_format",
        "{:.2f}".format,
    ):
        print(
            rules.drop(columns=["median_rss_mb", "median_min"]).to_string(index=False)
        )
    print(f"\nWrote {args.prefix}_jobs.tsv and {args.prefix}_rules.tsv")


if __name__ == "__main__":
    main()
