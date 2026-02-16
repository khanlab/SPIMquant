# Running Workflows

This guide covers running SPIMquant workflows in various environments.

## Workflow Overview

SPIMquant workflows consist of several stages:

1. **Data Loading**: Read BIDS-formatted SPIM data
2. **Registration**: Align to template using deformable registration
3. **Segmentation**: Detect pathology using various methods
4. **Quantification**: Extract statistics by brain region
5. **Visualization**: Generate quality control outputs

## Local Execution

### Single Subject

Process one subject at a time:

```bash
pixi run spimquant /bids /output participant \
  --filter_subjects 01 \
  --cores all
```

### Multiple Subjects

Process all subjects in parallel:

```bash
pixi run spimquant /bids /output participant --cores all
```

The `--cores` option controls parallelization:

- `--cores all`: Use all available cores
- `--cores N`: Use N cores
- `--cores 1`: Sequential processing

## Cluster Execution

<!-- TODO: Add detailed cluster execution guide -->

### SLURM

SPIMquant can submit jobs to SLURM clusters:

```bash
pixi run spimquant /bids /output participant \
  --profile slurm \
  --jobs 100 \
  --cores all
```

### Custom Profiles

<!-- TODO: Document how to create custom execution profiles -->

Create a custom profile in `spimquant/profiles/`:

```yaml
# Example SLURM profile
cluster: "sbatch -p {params.partition} -t {params.time}"
jobs: 100
```

## Cloud Execution

<!-- TODO: Expand cloud execution documentation -->

### Coiled

Run on cloud infrastructure with Coiled:

```bash
pixi run spimquant /bids /output participant \
  --cloud \
  --cores all
```

### Cloud Storage

Read data directly from cloud storage:

```bash
# S3
pixi run spimquant s3://bucket/bids /output participant --cores all

# GCS
pixi run spimquant gs://bucket/bids /output participant --cores all
```

## Workflow Control

### Dry Run

Test workflow without execution:

```bash
pixi run spimquant /bids /output participant -n
```

This shows:

- Jobs that would be executed
- Input/output files
- Rule dependencies
- Estimated resource usage

### Partial Execution

Run specific parts of the workflow:

```bash
# Run until specific rule
pixi run spimquant ... --until convert_to_nifti

# Run from specific rule
pixi run spimquant ... --from register_to_template

# Force re-run specific rule
pixi run spimquant ... --forcerun segment_pathology
```

### Resume After Failure

SPIMquant automatically resumes from the last successful step:

```bash
# Simply re-run the same command
pixi run spimquant /bids /output participant --cores all
```

## Monitoring Progress

### Real-time Monitoring

Watch progress during execution:

```bash
# Enable verbose output
pixi run spimquant ... --verbose

# Show executed commands
pixi run spimquant ... --printshellcmds
```

### Progress Tracking

<!-- TODO: Add progress tracking details -->

Snakemake displays:

- Completed jobs
- Running jobs
- Pending jobs
- Failed jobs

### Log Files

Check logs for detailed information:

```bash
# Snakemake log
cat /output/.snakemake/log/*.log

# Job-specific logs
cat /output/.snakemake/log/<job_name>/*.log
```

## Resource Management

### Memory Limits

Set memory constraints:

```bash
pixi run spimquant ... --resources mem_mb=64000
```

### Thread Control

Control threading per job:

```bash
pixi run spimquant ... --cores 16 --set-threads convert_to_nifti=4
```

### Time Limits

<!-- TODO: Add time limit configuration -->

Set maximum execution time for cluster jobs.

## Workflow Visualization

### DAG Visualization

Generate workflow graph:

```bash
pixi run spimquant /bids /output participant --dag | dot -Tpng > dag.png
```

### Rule Graph

Show rule dependencies:

```bash
pixi run spimquant /bids /output participant --rulegraph | dot -Tpng > rules.png
```

## Troubleshooting

### Failed Jobs

When jobs fail:

1. Check log files in `.snakemake/log/`
2. Re-run with `--verbose` for detailed output
3. Use `--printshellcmds` to see exact commands
4. Check input file existence and permissions

### Cleaning Up

Remove intermediate files:

```bash
# Clean specific outputs
pixi run spimquant ... --delete-temp-output

# Clean all outputs (use with caution)
pixi run spimquant ... --delete-all-output
```

### Lock Files

If workflow is interrupted, remove locks:

```bash
# Remove lock directory
rm -rf /output/.snakemake/locks/
```

## Best Practices

<!-- TODO: Add workflow best practices -->

1. **Always dry run first**: Use `-n` to validate
2. **Monitor resources**: Watch memory and CPU usage
3. **Use checkpoints**: Save intermediate results
4. **Clean temp files**: Free disk space after completion
5. **Generate reports**: Document your workflow execution

## Next Steps

- [Group Analysis](group_analysis.md): Statistical comparisons
- [Cloud Processing](cloud.md): Scale to cloud infrastructure
- [Examples](../examples/workflows.md): Complete workflow examples