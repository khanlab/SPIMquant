# Processing Cloud Data

<!-- TODO: Add cloud data processing guide -->

Process BIDS datasets stored in cloud storage.

## S3 Storage

<!-- TODO: Add S3 guide -->

```bash
pixi run spimquant s3://bucket/bids /output participant --cores all
```

## GCS Storage

<!-- TODO: Add GCS guide -->

```bash
pixi run spimquant gs://bucket/bids /output participant --cores all
```

## MIND

Datasets hosted on the MIND platform are read by `mind://` URL; see
[Processing Datasets on MIND](mind_data.md).

## Performance Tips

<!-- TODO: Add performance tips -->

## Next Steps

- [Cloud Processing](../usage/cloud.md)