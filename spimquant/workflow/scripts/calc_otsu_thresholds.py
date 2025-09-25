import json

import dask.array as da
import numpy as np
from skimage.filters import threshold_multiotsu

histogram = da.from_zarr(snakemake.params.histogram_uri).compute()

thresholds = dict()
for k in range(2, snakemake.params.otsu_max_k + 1):
    print(f"calculating otsu for k={k}")
    thresholds[k] = list()
    thresholds[k].append(0)
    thresholds[k].extend(threshold_multiotsu(hist=histogram, classes=k).tolist())
    thresholds[k].append(1e10)  # max
    print(thresholds[k])


with open(snakemake.output.otsu_thresholds, "w") as f:
    json.dump(thresholds, f, indent=4)
