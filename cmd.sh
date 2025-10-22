#!/bin/bash

export DASK_LOGGING=debug

source ./activate
spimquant /nfs/khan/datasets/MIND/mouse_app_lecanemab/ /tmp participant --participant-label AS36F4 -c all --use-conda --conda-prefix /home/UWO/akhan488/conda-snakemake-envs/ --segmentation-level 0 
cp -R /tmp/* /nfs/khan/datasets/MIND/mouse_app_lecanemab/derivatives/spimquant_v0.3.0rc
