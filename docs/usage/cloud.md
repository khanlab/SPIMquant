# Cloud Processing

SPIMquant supports cloud-based processing for scalable analysis of large datasets.

## Overview

Cloud processing enables:

- **Scalable compute**: Process many subjects in parallel
- **Cloud storage**: Direct access to S3, GCS data
- **Cost efficiency**: Pay only for resources used
- **No local hardware**: Process without local infrastructure

## Cloud Storage Support

### Amazon S3

Read BIDS datasets directly from S3:

```bash
pixi run spimquant s3://bucket-name/bids/ /local/output participant --cores all
```

#### S3 Configuration

Set AWS credentials:

```bash
export AWS_ACCESS_KEY_ID="your-access-key"
export AWS_SECRET_ACCESS_KEY="your-secret-key"
export AWS_DEFAULT_REGION="us-east-1"
```

Or use AWS CLI configuration:

```bash
aws configure
```

### Google Cloud Storage

Read from GCS buckets:

```bash
pixi run spimquant gs://bucket-name/bids/ /local/output participant --cores all
```

#### GCS Configuration

Authenticate with gcloud:

```bash
gcloud auth application-default login
```

Or use service account:

```bash
export GOOGLE_APPLICATION_CREDENTIALS="/path/to/service-account.json"
```

## Cloud Execution with Coiled

<!-- TODO: Add detailed Coiled setup and usage -->

SPIMquant integrates with [Coiled](https://coiled.io) for cloud execution.

### Setup Coiled

1. Create Coiled account at [coiled.io](https://coiled.io)
2. Install Coiled CLI:
   ```bash
   pip install coiled
   ```
3. Authenticate:
   ```bash
   coiled login
   ```

### Run on Cloud

```bash
pixi run spimquant /bids /output participant \
  --cloud \
  --cores all
```

### Coiled Configuration

<!-- TODO: Add Coiled configuration details -->

Configure cloud resources:

```yaml
# TODO: Add example Coiled config
```

## Hybrid Workflows

Combine local and cloud resources:

```bash
# Download from cloud, process locally
pixi run spimquant s3://bucket/bids /local/output participant --cores all

# Process locally, upload results to cloud
pixi run spimquant /local/bids s3://bucket/output participant --cores all
```

## Cost Optimization

<!-- TODO: Add cost optimization strategies -->

### Data Transfer Costs

Minimize data transfer:

- Process in same region as data
- Use cloud storage classes appropriately
- Clean up intermediate files

### Compute Costs

Optimize compute resources:

- Right-size instance types
- Use spot instances when possible
- Stop resources when not in use

### Storage Costs

Manage storage efficiently:

- Delete temporary files after completion
- Archive old results to cheaper storage tiers
- Use compression for large outputs

## Cloud Provider Guides

### AWS

<!-- TODO: Add AWS-specific guide -->

#### Setting Up

1. Create S3 bucket
2. Configure IAM permissions
3. Set up EC2 instances (if not using Coiled)

#### Running on AWS

```bash
# TODO: Add AWS execution examples
```

### Google Cloud Platform

<!-- TODO: Add GCP-specific guide -->

#### Setting Up

1. Create GCS bucket
2. Configure service account
3. Set up Compute Engine instances

#### Running on GCP

```bash
# TODO: Add GCP execution examples
```

### Azure

<!-- TODO: Add Azure support documentation -->

Azure support is planned for future releases.

## Monitoring and Debugging

### Cloud Logs

<!-- TODO: Add cloud logging details -->

Access logs for cloud runs:

```bash
# TODO: Add log access commands
```

### Resource Monitoring

Monitor cloud resource usage:

- CPU utilization
- Memory usage
- Network I/O
- Storage I/O

### Troubleshooting

Common cloud issues:

1. **Authentication failures**: Check credentials
2. **Permission errors**: Verify IAM/service account permissions
3. **Region errors**: Ensure resources in same region
4. **Network timeouts**: Increase timeout settings

## Security Considerations

<!-- TODO: Add security best practices -->

### Data Security

- Use encrypted storage buckets
- Enable encryption in transit
- Implement access controls
- Audit access logs

### Credential Management

- Never commit credentials to code
- Use environment variables or secret management
- Rotate credentials regularly
- Use least-privilege access

## Performance Tips

### Network Performance

- Process data in same region/zone
- Use high-bandwidth instances
- Enable accelerated networking

### Storage Performance

- Use SSD-backed storage
- Enable caching where appropriate
- Parallelize I/O operations

### Compute Performance

- Choose appropriate instance types
- Use instances with local SSD for temp files
- Enable instance-level parallelization

## Example Workflows

### Complete Cloud Workflow

```bash
# 1. Upload data to S3
aws s3 sync /local/bids s3://bucket/bids/

# 2. Process on cloud
pixi run spimquant s3://bucket/bids s3://bucket/output participant \
  --cloud \
  --cores all

# 3. Download results
aws s3 sync s3://bucket/output /local/output/
```

### Hybrid Processing

```bash
# Process participant-level on cloud
pixi run spimquant s3://bucket/bids /local/output participant --cloud --cores all

# Run group-level locally
pixi run spimquant /local/bids /local/output group \
  --contrast_column treatment \
  --contrast_values control drug
```

## Comparison: Local vs Cloud

| Aspect | Local | Cloud |
|--------|-------|-------|
| Setup | Hardware required | Account setup only |
| Cost | Upfront hardware | Pay-as-you-go |
| Scalability | Limited by hardware | Unlimited scaling |
| Data transfer | None | Can be significant |
| Maintenance | Manual | Managed |

## Next Steps

- [Configuration](configuration.md): Configure cloud storage
- [Workflows](workflows.md): Execution strategies
- [Examples](../examples/workflows.md): Complete cloud examples