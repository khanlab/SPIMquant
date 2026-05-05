"""Create a BIDS JSON sidecar file that describes the columns of a tabular output.

This script is called by the generic ``create_tabular_json_sidecar`` Snakemake rule.
It reads the column names from a TSV or Parquet file, looks each one up in the
``column_descriptions`` dictionary supplied via ``snakemake.params``, and writes a
BIDS-compliant JSON sidecar next to the tabular file.

Column name resolution (applied in order):
1. Exact match in ``column_descriptions``.
2. Stain-prefixed compound column: ``{stain}+{metric}`` — the metric description is
   looked up and the LongName / Description are annotated to mention the stain.
3. Stat-suffixed column: ``{base}_{stat}`` where *stat* is one of the known
   ``stats_map`` suffixes (``tstat``, ``pval``, ``cohensd``) — description is derived
   from the base column and the stat entry.
4. Group-stat column with group label: ``{base}_{stat}_{group}`` — same as above but
   also mentions the group label.
5. Subject-count column: ``n_{group}`` — generic description referencing the group.
6. Mean/std/min/max aggregate columns: ``{base}_mean``, ``{base}_std``, etc.
7. If no match is found, the column is included without a description so that the JSON
   remains complete (all columns are listed).

This is a Snakemake script that expects the ``snakemake`` object to be available.
"""

import json
import os
import re
from pathlib import Path


def get_columns(filepath):
    """Return list of column names from a TSV or Parquet file."""
    ext = Path(filepath).suffix.lower()
    if ext == ".parquet":
        import pyarrow.parquet as pq

        schema = pq.read_schema(filepath)
        return list(schema.names)
    else:
        # TSV – read only the header row for speed
        with open(filepath) as fh:
            header = fh.readline().rstrip("\n")
        return header.split("\t")


def resolve_column(col, column_descriptions, stat_suffixes):
    """Return a BIDS-style description dict for *col*, or an empty dict if unknown.

    Parameters
    ----------
    col:
        Column name to resolve.
    column_descriptions:
        Mapping of known column names to BIDS description dicts (from config).
    stat_suffixes:
        Set of recognised statistical suffix tokens (e.g. ``{'tstat', 'pval', ...}``).
    """
    # 1. Exact match
    if col in column_descriptions:
        return dict(column_descriptions[col])

    # 2. Stain-prefixed compound: "{stain}+{metric}"
    if "+" in col:
        stain, _, metric = col.partition("+")
        # Allow nested stat suffixes in the metric part (e.g. "Abeta+fieldfrac_tstat")
        metric_desc = resolve_column(metric, column_descriptions, stat_suffixes)
        if metric_desc:
            entry = dict(metric_desc)
            long_name = entry.get("LongName", metric)
            description = entry.get("Description", "")
            entry["LongName"] = f"{long_name} ({stain})"
            entry["Description"] = (
                f"{description} Stain/channel: {stain}."
                if description
                else f"Metric '{metric}' for stain/channel '{stain}'."
            )
            return entry
        return {
            "LongName": f"{metric} ({stain})",
            "Description": f"Metric '{metric}' for stain/channel '{stain}'.",
        }

    # 3 & 4. Stat-suffixed column: "{base}_{stat}" or "{base}_{stat}_{group}"
    #   Walk through known stat suffixes to find a match.
    for stat in stat_suffixes:
        # Pattern: base_stat_group  (group may contain underscores)
        pattern_with_group = re.compile(rf"^(.+)_{re.escape(stat)}_(.+)$")
        m = pattern_with_group.match(col)
        if m:
            base, group = m.group(1), m.group(2)
            base_desc = resolve_column(base, column_descriptions, stat_suffixes)
            stat_desc = column_descriptions.get(stat, {})
            long_name = base_desc.get("LongName", base) if base_desc else base
            stat_long = stat_desc.get("LongName", stat)
            return {
                "LongName": f"{long_name} – {stat_long} (group: {group})",
                "Description": (
                    f"{stat_desc.get('Description', stat_long)} "
                    f"for metric '{base}', group '{group}'."
                ),
            }
        # Pattern: base_stat
        pattern_no_group = re.compile(rf"^(.+)_{re.escape(stat)}$")
        m = pattern_no_group.match(col)
        if m:
            base = m.group(1)
            base_desc = resolve_column(base, column_descriptions, stat_suffixes)
            stat_desc = column_descriptions.get(stat, {})
            long_name = base_desc.get("LongName", base) if base_desc else base
            stat_long = stat_desc.get("LongName", stat)
            return {
                "LongName": f"{long_name} – {stat_long}",
                "Description": (
                    f"{stat_desc.get('Description', stat_long)} for metric '{base}'."
                ),
            }

    # 5. Aggregate suffix: "{base}_mean", "{base}_std", "{base}_min", "{base}_max"
    agg_suffixes = ("_mean", "_std", "_min", "_max")
    for agg in agg_suffixes:
        if col.endswith(agg):
            base = col[: -len(agg)]
            base_desc = resolve_column(base, column_descriptions, stat_suffixes)
            agg_label = agg.lstrip("_").capitalize()
            if base_desc:
                long_name = base_desc.get("LongName", base)
                desc_text = base_desc.get("Description", "")
                entry = {
                    "LongName": f"{long_name} – {agg_label}",
                    "Description": f"{agg_label} of: {desc_text}"
                    if desc_text
                    else f"{agg_label} of '{base}'.",
                }
                if "Units" in base_desc:
                    entry["Units"] = base_desc["Units"]
                return entry
            return {
                "LongName": f"{base} – {agg_label}",
                "Description": f"{agg_label} of '{base}'.",
            }

    # 6. Subject-count column: "n_{group}"
    if re.match(r"^n_[A-Za-z0-9_]+$", col):
        group = col[2:]
        if group == "subjects":
            return dict(column_descriptions.get("n_subjects", {
                "LongName": "Number of subjects",
                "Description": "Total number of subjects contributing data to this row.",
            }))
        return {
            "LongName": f"Number of subjects (group: {group})",
            "Description": (
                f"Number of subjects in group '{group}' contributing data to this row."
            ),
        }

    # No match – return empty dict; column will still appear in JSON without metadata
    return {}


def build_sidecar(columns, column_descriptions, stat_suffixes):
    """Build the full BIDS sidecar dict for all *columns*."""
    sidecar = {}
    for col in columns:
        entry = resolve_column(col, column_descriptions, stat_suffixes)
        sidecar[col] = entry
    return sidecar


def main():
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    column_descriptions = snakemake.params.column_descriptions
    # stats_maps contains the short stat tokens used as suffixes (tstat, pval, cohensd)
    stat_suffixes = set(snakemake.params.get("stats_maps", ["tstat", "pval", "cohensd"]))

    columns = get_columns(input_file)
    sidecar = build_sidecar(columns, column_descriptions, stat_suffixes)

    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(output_file, "w") as fh:
        json.dump(sidecar, fh, indent=4)


main()
