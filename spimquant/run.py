#!/usr/bin/env python3
from pathlib import Path

from snakebids import bidsapp, plugins
import os
import shutil


app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.CliConfig(),
        plugins.BidsArgs(bids_dir=False),
#        plugins.BidsValidator(),
        plugins.Version(distribution="spimquant"),
    ]
)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
