# Hardware Requirements

<!-- TODO: This section needs careful editing - some specific constraints and benchmarks may need verification -->

SPIMquant processes large-scale microscopy datasets and requires substantial computational resources. This guide helps you plan your hardware setup.

## Minimum Requirements

For processing typical SPIM datasets:

- **CPU**: 8+ cores
- **RAM**: 32 GB
- **Storage**: 500 GB free space (SSD recommended)
- **OS**: Linux (64-bit)

## Recommended Configuration

For optimal performance:

- **CPU**: 16-32 cores
- **RAM**: 64-128 GB
- **Storage**: 1+ TB fast SSD
- **OS**: Linux (Ubuntu 20.04+ or similar)

## Resource Considerations

### Memory (RAM)

Memory requirements depend on:

1. **Dataset Resolution**: Higher resolution requires more memory
2. **Downsampling Level**: Setting a higher downsampling level for registration or segmentation will reduce memory requirements significantly
3. **Parallel Processing**: More CPU's available to Dask = more memory needed, can use --set-threads to customize per rule

**Memory Guidelines:**

| Task | Memory per Core | Example |
|------|----------------|---------|
| Low-res processing | 2-4 GB | 8 cores = 16-32 GB |
| High-res processing | 4-8 GB | 8 cores = 32-64 GB |
| Template registration | 16-32 GB | Single job |


### CPU Cores

SPIMquant parallelizes at multiple levels:

1. **Snakemake Level**: Multiple subjects/samples processed concurrently
2. **Dask Level**: Within-job parallelization for array operations

**Core Guidelines:**

- **Minimum**: 4 cores
- **Recommended**: 16+ cores
- **High-throughput**: 32+ cores

More cores significantly reduce processing time for multi-subject datasets.

### Storage

Storage needs vary by dataset size:

| Component | Space Needed | Notes |
|-----------|-------------|-------|
| Input data | 1-2000 GB/subject | Depends on resolution |
| Intermediate files | 2-3x input | Temporary images |
| Output results | 0.1x input | Registered data + stats |

**Storage Recommendations:**

- **SSD**: Strongly recommended for temp files
- **Network Storage**: Possible but slower
- **Cleanup**: Intermediate files are automatically deleted

### Disk I/O

Fast storage improves performance:

- **SSD**: 5-10x faster than HDD
- **NVMe**: Best performance for intensive workflows
- **Network**: Can bottleneck if many chunks

Consider local storage for:

- Temporary files
- Frequently accessed data
- Registration intermediate steps

## Scaling Strategies

### Single Workstation

For local processing:

```bash
# Use all cores, but control memory
pixi run spimquant ... --cores all --resources mem_mb=64000
```

### HPC Cluster

For SLURM or other clusters:

```bash
# Submit jobs to cluster each job using 32 cores, 100 jobs submitted at a time
pixi run spimquant ... \
  --profile cluster \
  --jobs 100 \
  --cores 32
```

See [Cluster Configuration](../usage/workflows.md#cluster-execution) for details.

### Cloud Processing

For cloud-based processing see [Cloud Processing Guide](../usage/cloud.md) for setup.

## Performance Tuning


### Monitor Resource Usage

During processing, monitor resources:

```bash
# CPU and memory usage
htop

# Disk I/O
iotop

# Disk space
df -h
```

### Memory Management Tips

1. **Process in batches**: Run subset of subjects if memory-limited, can use --particiant-label or --filter-spim
2. **Reduce resolution**: Use lower-resolution pyramids for initial testing (--registration-level, --segmentation-level)
3. **Limit concurrent jobs**: Use `--cores` to control parallelism
4. **Close other applications**: Free up memory for SPIMquant
5. **Skip segmentation if you only need atlas registration**: use --no-segmentation flag


## Benchmarks

Typical processing times on different configurations:

(placeholder)

| Configuration | Dataset | Time |
|--------------|---------|------|
| 8 cores, 32 GB | 1 subject |  _ hours |
| 16 cores, 64 GB | 5 subjects  | _ hours |
| 32 cores, 128 GB | 20 subjects | _ hours |

!!! note "Benchmarks"
    Times vary significantly based on dataset resolution, template choice, and specific workflow steps enabled.


## Getting Help

If you're unsure about hardware requirements for your specific use case:

1. Start with minimum requirements
2. Monitor resource usage during test runs
3. Scale up as needed
4. Ask on [GitHub Discussions](https://github.com/khanlab/SPIMquant/discussions)