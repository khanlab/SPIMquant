"""
Run Atropos segmentation with automatic k-retry logic.

Tries GMM segmentation with init_k classes, and decrements k if Atropos fails,
down to min_k. This handles cases where the image doesn't have enough distinct
intensity classes to support the initial k value.
"""

import os
import shutil
import subprocess

init_k = snakemake.params.init_k
min_k = snakemake.params.min_k
threads = snakemake.threads
mrf_smoothing = snakemake.params.mrf_smoothing
mrf_radius = snakemake.params.mrf_radius

downsampled = snakemake.input.downsampled
mask = snakemake.input.mask
dseg = snakemake.output.dseg
posteriors_dir = snakemake.output.posteriors_dir

os.makedirs(posteriors_dir, exist_ok=True)

for k in range(init_k, min_k - 1, -1):
    cmd = (
        f"ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} "
        f"Atropos -v -d 3 --initialization KMeans[{k}] "
        f" --intensity-image {downsampled} "
        f" --output [{dseg},{posteriors_dir}/class-%02d.nii] "
        f" --mask-image {mask} --mrf [{mrf_smoothing},{mrf_radius}]"
    )
    result = subprocess.run(cmd, shell=True)
    if result.returncode == 0:
        break
    elif k == min_k:
        raise subprocess.CalledProcessError(result.returncode, cmd)
    else:
        # Clean up and retry with lower k
        if os.path.exists(dseg):
            os.remove(dseg)
        shutil.rmtree(posteriors_dir, ignore_errors=True)
        os.makedirs(posteriors_dir, exist_ok=True)
