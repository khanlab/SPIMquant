# SPIMquant
[![Documentation Status](https://readthedocs.org/projects/spimquant/badge/?version=latest)](https://spimquant.readthedocs.io/en/latest/?badge=latest)

SPIMquant takes lightsheet (SPIM) scans of whole cleared brains, organized as
a BIDS dataset, registers them to a template, segments each stain at full
resolution with a method chosen per stain, and summarises the results per
atlas region, per object and across a cohort, with QC reports at every step.
It is a Snakemake workflow packaged as a BIDS App, and runs on a workstation
or a SLURM cluster.

**Documentation: <https://spimquant.readthedocs.io>**, including
[example outputs](https://spimquant.readthedocs.io/en/latest/examples/) and
the [command-line reference](https://spimquant.readthedocs.io/en/latest/reference/cli/).

## Quick start

```bash
curl -fsSL https://pixi.sh/install.sh | bash      # if pixi is not installed
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant
pixi run spimquant /path/to/bids /path/to/output participant -n   # dry run
pixi run spimquant /path/to/bids /path/to/output participant --cores all
pixi run spimquant /path/to/bids /path/to/output group --cores all
```

Use `pixi run -e gpu ...` on a machine with an NVIDIA GPU for deep-learning
segmentation (LANTERN plaques, vesselFM vessels). SPIMquant is also published
as a conda package; see
[Installation](https://spimquant.readthedocs.io/en/latest/getting_started/installation/).

## Contributing

See [Contributing](https://spimquant.readthedocs.io/en/latest/contributing/)
and, for coding agents, [AGENTS.md](AGENTS.md). Checks:

```bash
pixi run -e dev quality_check
pixi run -e dev pytest tests/
```

## License

MIT. See [LICENSE](LICENSE).
