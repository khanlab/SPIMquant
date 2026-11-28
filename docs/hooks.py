"""MkDocs hooks: generate the reference pages before each build.

See docs/gen/build.py for what is generated and from where, and AGENTS.md
for which parts of the docs are maintained by hand.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "gen"))

import build  # noqa: E402


def on_pre_build(config):
    # Fail the build rather than publish stale or missing reference pages.
    build.generate()
