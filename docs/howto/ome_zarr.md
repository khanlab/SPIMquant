# Working with OME-Zarr

<!-- TODO: Add OME-Zarr guide -->

This guide covers working with OME-Zarr files in SPIMquant.

## OME-Zarr Format

OME-Zarr is a cloud-optimized format for storing multidimensional microscopy data.

## Reading OME-Zarr Data

<!-- TODO: Add reading instructions -->

```bash
pixi run spimquant /bids /output participant \
  --filter-spim extension='ome.zarr.zip' \
  --cores all
```

## Zarr Zipstores

<!-- TODO: Add zipstore information -->

## Creating OME-Zarr Files

<!-- TODO: Add creation guide -->

See [zarrnii](https://github.com/khanlab/zarrnii) for creating OME-Zarr files.

## Next Steps

- [Segmentation Methods](segmentation.md)
- [Cloud Data](cloud_data.md)