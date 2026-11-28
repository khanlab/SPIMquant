#!/usr/bin/env python3
"""Regenerate docs/generated/ from dry runs of the workflow on tests/bids_ds.

Run with ``pixi run docs-introspect`` (needs the workflow environment). For
each scenario in docs/gen/scenarios.yml it builds the job DAG, without running
anything, and records:

- ``docs/generated/rulegraph_<scenario>.dot``: the rule graph, with nodes and
  edges sorted so that the file only changes when the workflow does;
- ``docs/generated/outputs.tsv``: every file a job writes, with the entity
  values of the test dataset replaced by placeholders (``sub-{subject}``,
  ``stain-{stain}``...), the rule that writes it, whether it is temporary, and
  the first scenario that produces it.

CI runs this and fails if the committed files differ, so they cannot drift
from the workflow. docs/gen/build.py turns them into the reference pages.
"""

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
GEN = REPO / "docs" / "generated"
SCENARIOS = Path(__file__).with_name("scenarios.yml")

# --- worker: runs inside one spimquant process ------------------------------


def _worker(out_json, argv):
    """Patch --rulegraph to dump the DAG as JSON, then run spimquant."""
    import snakemake.api  # noqa: F401  (import order avoids a circular import)
    from snakemake.io import is_flagged
    from snakemake.workflow import Workflow

    def dump(self):
        self._prepare_dag(forceall=False, ignore_incomplete=True, lock_warn_only=True)
        self._build_dag()
        jobs = [
            {
                "rule": job.rule.name,
                "outputs": [[str(f), bool(is_flagged(f, "temp"))] for f in job.output],
            }
            for job in self.dag.jobs
        ]
        Path(out_json).write_text(
            json.dumps({"rulegraph": self.dag.rule_dot(), "jobs": jobs})
        )

    Workflow.printrulegraph = dump
    sys.path.insert(0, str(REPO))
    from spimquant.run import app

    sys.argv = ["spimquant", *argv, "--rulegraph"]
    app.run()


# --- normalisation ----------------------------------------------------------


def normalize_dot(dot):
    """Sort node and edge lines so the graph text is deterministic."""
    lines = dot.strip().splitlines()
    head = [l for l in lines if "->" not in l and not re.match(r"\s*\d+\[", l)]
    nodes = sorted(l.strip() for l in lines if re.match(r"\s*\d+\[", l))
    edges = sorted({l.strip() for l in lines if "->" in l})
    # head ends with the closing brace
    close = head.pop() if head and head[-1].strip() == "}" else "}"
    return (
        "\n".join(
            head + ["\t" + n for n in nodes] + ["\t" + e for e in edges] + [close]
        )
        + "\n"
    )


def make_normalizer(entities):
    """Build a function that replaces test-dataset entity values with
    placeholders, from the ``placeholders`` mapping in scenarios.yml."""
    subs = []
    for name, values in entities.items():
        alt = "|".join(re.escape(v) for v in sorted(values, key=len, reverse=True))
        if name == "stain":
            subs.append((re.compile(rf"(stain-)({alt})(?=[_.])"), r"\1{stain}"))
            subs.append((re.compile(rf"(?<=[_/])({alt})\+"), "{stain}+"))
        elif name == "template":
            subs.append(
                (
                    re.compile(rf"(?<![^_/])(space|from|to|tpl)-({alt})(?=[_./])"),
                    r"\1-{template}",
                )
            )
        else:
            entity = {
                "subject": "sub",
                "atlas": "seg",
                "seg_method": "desc",
                "object_method": "obj",
            }.get(name, name)
            subs.append(
                (
                    re.compile(rf"(?<![^_/]){entity}-({alt})(?=[_./])"),
                    f"{entity}-{{{name}}}",
                )
            )

    def norm(path):
        path = re.sub(r"^group/[^/]+/", "group/{label}/", path)
        path = re.sub(r"(?<=_)metric-[^_]+(?=_)", "metric-{metric}", path)
        for rx, rep in subs:
            path = rx.sub(rep, path)
        return path

    return norm


# --- main -------------------------------------------------------------------


def run_scenario(sc, workdir):
    bids = REPO / sc.get("bids", "tests/bids_ds")
    out = Path(workdir) / sc["name"]
    out_json = Path(workdir) / f"{sc['name']}.json"
    prep = sc.get("prepare")
    if prep:  # e.g. add a second subject for the cohort template
        ds = Path(workdir) / f"bids_{sc['name']}"
        subprocess.run(["cp", "-r", str(bids), str(ds)], check=True)
        subprocess.run(["bash", "-c", prep], cwd=ds, check=True)
        bids = ds
    argv = [str(bids), str(out), sc.get("level", "participant"), *sc.get("args", [])]
    proc = subprocess.run(
        [sys.executable, __file__, "--worker", str(out_json), *argv],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not out_json.exists():
        sys.exit(
            f"scenario {sc['name']} failed:\n{proc.stdout[-3000:]}\n{proc.stderr[-3000:]}"
        )
    return json.loads(out_json.read_text())


def main():
    cfg = yaml.safe_load(SCENARIOS.read_text())
    norm = make_normalizer(cfg["placeholders"])
    GEN.mkdir(parents=True, exist_ok=True)
    rows = {}
    with tempfile.TemporaryDirectory() as tmp:
        for sc in cfg["scenarios"]:
            print(f"scenario {sc['name']}: {' '.join(sc.get('args', []))}", flush=True)
            d = run_scenario(sc, tmp)
            if sc.get("rulegraph", True):
                (GEN / f"rulegraph_{sc['name']}.dot").write_text(
                    normalize_dot(d["rulegraph"])
                )
            for job in d["jobs"]:
                for f, temp in job["outputs"]:
                    if f.startswith(("/", "..")):  # outside the output dir
                        continue
                    key = norm(f)
                    if key not in rows:
                        rows[key] = (job["rule"], temp, sc["name"])
    lines = ["path\trule\ttemp\tscenario"]
    for path, (rule, temp, scen) in sorted(rows.items()):
        lines.append(f"{path}\t{rule}\t{'yes' if temp else 'no'}\t{scen}")
    (GEN / "outputs.tsv").write_text("\n".join(lines) + "\n")
    print(f"wrote {len(rows)} output patterns to {GEN / 'outputs.tsv'}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--worker":
        _worker(sys.argv[2], sys.argv[3:])
    else:
        main()
