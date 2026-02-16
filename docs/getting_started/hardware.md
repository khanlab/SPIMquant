# Hardware Requirements

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
2. **Registration Method**: `greedy` deformable registration is memory-intensive
3. **Parallel Processing**: More concurrent jobs = more memory needed

**Memory Guidelines:**

| Task | Memory per Core | Example |
|------|----------------|---------|
| Low-res processing | 2-4 GB | 8 cores = 16-32 GB |
| High-res processing | 4-8 GB | 8 cores = 32-64 GB |
| Template registration | 16-32 GB | Single job |

!!! warning "Out of Memory"
    The `greedy` registration step can require 16-32 GB for a single job. Ensure sufficient memory is available or reduce concurrent jobs.

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
| Input data | 1-100 GB/subject | Depends on resolution |
| Intermediate files | 2-5x input | Temporary registrations |
| Output results | 0.5-2x input | Registered data + stats |

**Storage Recommendations:**

- **SSD**: Strongly recommended for temp files
- **Network Storage**: Possible but slower
- **Cleanup**: Intermediate files can be deleted

!!! tip "Storage Management"
    Use `--delete-temp-output` to automatically clean intermediate files after successful completion.

### Disk I/O

Fast storage improves performance:

- **SSD**: 5-10x faster than HDD
- **NVMe**: Best performance for intensive workflows
- **Network**: Can bottleneck on large files

Consider local storage for:

- Temporary files
- Frequently accessed data
- Registration intermediate steps

## Scaling Strategies

### Single Workstation

For local processing:

```bash
# Use all cores, but control memory
pixi run spimquant ... --cores 16 --resources mem_mb=64000
```

### HPC Cluster

For SLURM or other clusters:

```bash
# Submit jobs to cluster
pixi run spimquant ... \
  --profile cluster \
  --jobs 100 \
  --cores all
```

See [Cluster Configuration](../usage/workflows.md#cluster-execution) for details.

### Cloud Processing

For cloud-based processing:

```bash
# Use Coiled for cloud execution
pixi run spimquant ... \
  --cloud \
  --cores all
```

See [Cloud Processing Guide](../usage/cloud.md) for setup.

## Performance Tuning

### Optimize for Your Hardware

**High Memory, Few Cores:**
```bash
# Run fewer jobs in parallel, but use more threads per job
pixi run spimquant ... --cores 4 --threads 8
```

**Many Cores, Limited Memory:**
```bash
# Run many light jobs in parallel
pixi run spimquant ... --cores all --threads 2
```

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

1. **Process in batches**: Run subset of subjects if memory-limited
2. **Reduce resolution**: Use lower-resolution pyramids for initial testing
3. **Limit concurrent jobs**: Use `--cores` to control parallelism
4. **Close other applications**: Free up memory for SPIMquant

## Benchmarks

Typical processing times on different configurations:

| Configuration | Dataset | Time |
|--------------|---------|------|
| 8 cores, 32 GB | 1 subject, 10 GB | 2-4 hours |
| 16 cores, 64 GB | 5 subjects, 50 GB | 6-10 hours |
| 32 cores, 128 GB | 20 subjects, 200 GB | 12-24 hours |

!!! note "Benchmarks"
    Times vary significantly based on dataset resolution, template choice, and specific workflow steps enabled.

## Cloud vs. Local Processing

### Local Processing

**Pros:**
- No data transfer costs
- Full control over resources
- No cloud setup required

**Cons:**
- Limited by local hardware
- Requires significant upfront hardware investment
- Manual resource management

### Cloud Processing

**Pros:**
- Scalable resources on-demand
- No hardware maintenance
- Pay only for what you use

**Cons:**
- Data transfer time and costs
- Cloud setup complexity
- Ongoing operational costs

See [Cloud Processing Guide](../usage/cloud.md) for cloud deployment details.

## Getting Help

If you're unsure about hardware requirements for your specific use case:

1. Start with minimum requirements
2. Monitor resource usage during test runs
3. Scale up as needed
4. Ask on [GitHub Discussions](https://github.com/khanlab/SPIMquant/discussions)