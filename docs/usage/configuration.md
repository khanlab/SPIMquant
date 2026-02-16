# Configuration

SPIMquant is configured through a combination of command-line arguments and a YAML configuration file.

## Configuration File Location

The main configuration file is located at:

```
spimquant/config/snakebids.yml
```

## Basic Configuration

<!-- TODO: Add detailed configuration options -->

The configuration file defines:

- Input data specifications
- Analysis parameters
- Template settings
- Segmentation methods
- Output options

## Configuration Structure

```yaml
# TODO: Add example configuration structure

# BIDS input specification
pybids_inputs:
  spim:
    filters:
      suffix: 'SPIM'
      extension: '.ome.zarr'
    wildcards:
      - subject
      - sample
      - stain

# Analysis parameters
registration:
  template: 'ABAv3'
  method: 'greedy'
  
segmentation:
  method: 'otsu'
  stains:
    - 'Abeta'
    - 'Iba1'
```

## Template Configuration

<!-- TODO: Document template configuration options -->

SPIMquant supports multiple templates:

- **ABAv3**: Allen Brain Atlas v3
- **gubra**: Gubra atlas
- **MBMv3**: Marmoset Brain Maps v3
- **turone**: Turone atlas
- **MouseIn**: Mouse MRI template

## Segmentation Configuration

<!-- TODO: Document segmentation configuration -->

Configure segmentation methods and parameters:

```yaml
# Example configuration
segmentation:
  method: 'threshold'  # or 'otsu+k3i2'
  intensity_correction: 'gaussian'  # or 'n4'
```

## Advanced Configuration

### Dask Configuration

<!-- TODO: Add Dask configuration details -->

Configure parallel processing with Dask:

```yaml
dask:
  scheduler: 'threads'
  num_workers: 8
```

### Cloud Storage Configuration

<!-- TODO: Add cloud storage configuration -->

Configure S3 or GCS access:

```yaml
storage:
  s3:
    endpoint: 's3.amazonaws.com'
  gcs:
    project: 'my-project'
```

## Example Configurations

### Basic Processing

<!-- TODO: Add example configurations -->

```yaml
# Minimal configuration for basic processing
```

### Advanced Workflow

```yaml
# Configuration for advanced processing with custom parameters
```

## Configuration Validation

SPIMquant validates configuration on startup. Common validation errors:

<!-- TODO: Add common validation errors and fixes -->

- Missing required fields
- Invalid template names
- Incorrect file paths

## Next Steps

- [Running Workflows](workflows.md): Use your configuration
- [CLI Reference](cli.md): Override config with command-line options
- [Examples](../examples/workflows.md): See complete configurations