# Installation

SPIMquant runs on 64-bit Linux. Its dependencies, including greedy, ANTs,
c3d and PyTorch, are conda packages managed by [pixi](https://pixi.sh).

## From a clone of the repository

Use this to run the test dataset, the development tools, or the GPU
environment.

```bash
curl -fsSL https://pixi.sh/install.sh | bash   # if pixi is not installed
git clone https://github.com/khanlab/SPIMquant.git
cd SPIMquant
pixi run spimquant --help                       # installs the environment on first use
```

The repository defines several environments; choose one with `pixi run -e <env>`:

| Environment | For |
|---|---|
| `default` | Running SPIMquant on CPU. |
| `gpu` | Running with CUDA PyTorch, for LANTERN (`--seg_method lantern`), vesselFM and FireANTs. Needs an NVIDIA driver; on a GPU-less SLURM login node see [Running SPIMquant](../usage/running.md#gpu-steps). |
| `dev` | `default` plus formatters, pytest and napari. |
| `docs` | Building this documentation. |

## As a conda package

To only run SPIMquant, install the package from the three channels it is
published on:

```bash
pixi global install \
  --channel https://prefix.dev/conda-forge \
  --channel https://prefix.dev/bioconda \
  --channel https://prefix.dev/akhanf/khanlab \
  spimquant
spimquant --help
```

conda, mamba and micromamba work with the same channels.

## Disk and memory

For the example brain (339 GB of Imaris at 1.625 × 1.625 × 2.75 µm, Abeta
with LANTERN and CD31 with vesselFM) the output directory takes 110 GB. Most
of it (101 GB) is the full-resolution distance map to the vessels
(`vessels/*_dist.ozx`); without vessels (`--no_vessels`) it is under 10 GB.
Temporary files written during the run come on top of this; put them on fast
local disk with `--work_dir`. Memory per job is set by each rule from the
image size; [Running SPIMquant](../usage/running.md#resources) lists measured
run times and memory.

## Next

[Quick start](quickstart.md): run the test dataset.
