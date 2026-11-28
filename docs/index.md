# SPIMquant

SPIMquant takes lightsheet (SPIM) scans of whole cleared brains, organized as
a [BIDS](https://bids-specification.readthedocs.io/) dataset, and produces
atlas-based measurements of what each stain shows: amyloid plaques,
microglia, vessels, neurons, fibres. It is a [Snakemake](https://snakemake.readthedocs.io)
workflow packaged as a BIDS App, so one command processes one brain or a
cohort, on a workstation or a SLURM cluster.

For each brain it:

1. registers the scan to a template (the Allen CCFv3 by default, or a
   lightsheet template built from your own cohort);
2. segments each stain at full resolution, with a method chosen per stain
   (thresholds, a Gaussian-mixture tail, difference-of-Gaussians, or deep
   learning: LANTERN for plaques, vesselFM for vessels);
3. splits the masks into objects and measures them: size, intensity,
   position in subject and template space, distance to the nearest object of
   every other stain and to the nearest vessel, colocalization;
4. summarises everything per atlas region (field fraction, count, density,
   object properties) in one table per brain, and paints the values into the
   template;
5. writes QC reports for registration and segmentation.

At the group level it stacks the tables of all brains with the columns of
`participants.tsv`, fits a linear model per region and metric, and writes
contrast tables, heatmaps and statistical maps.

[See what the outputs look like](examples/index.md){ .md-button }
[Install](getting_started/installation.md){ .md-button }

## Where to start

| I want to... | Read |
|---|---|
| Install and run the test dataset | [Installation](getting_started/installation.md), [Quick start](getting_started/quickstart.md) |
| Prepare my own data | [Input data](getting_started/input_data.md) |
| Understand the processing | [How it works](workflow_overview.md), [workflow graphs](reference/workflow.md) |
| Choose segmentation methods | [Segmentation](howto/segmentation.md) |
| Run on a cluster or a GPU node | [Running SPIMquant](usage/running.md) |
| Compare groups | [Group analysis](usage/group_analysis.md) |
| Look up an option or a file | [Command line](reference/cli.md), [output files](reference/outputs.md) |

## Citing

A paper describing SPIMquant is in preparation. Until then, cite the
repository: <https://github.com/khanlab/SPIMquant>.

SPIMquant is released under the MIT License.
