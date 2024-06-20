"""
This file handles file persistence of common objects used in pipeline. For example, dataset references
that specifies the files in a dataset can be stored as a file and read using the interfaces in this
file. Note that many file format is a directory containing a single file, which may seem that the directory
is unnecessary but this is to accommodate for the possible need to add more files to store information
about the object class' instance.
"""


import shutil
import json
from .fs import ensure_dir_exists
from .strenc import get_encoder, get_decoder_hook
from .dataset_reference import DatasetReference


def write_dataset_reference(ref: DatasetReference, path: str):
    ensure_dir_exists(path, True)
    with open(f'{path}/dataset.json', 'w') as outfile:
        json.dump(ref, outfile, cls=get_encoder(), indent=2)


def read_dataset_reference(path: str) -> DatasetReference:
    with open(f'{path}/dataset.json', 'r') as infile:
        ref = json.load(infile, object_hook=get_decoder_hook())
    return ref


def write_cellseg3d_config(config, path: str, write_model=False):
    ensure_dir_exists(path, True)
    if write_model:
        cur_path = config.model_weight_path
        if cur_path is None:
            print('WARNING: cur_path is None while write_model is True!')
            config.model_weight_path = f'{path}/wnet.pth'
            shutil.copy(cur_path, config.model_weight_path)
    with open(f'{path}/model_config.json', 'w') as outfile:
        json.dump(ref, outfile, cls=get_encoder(), indent=2)


def read_cellseg3d_config(path: str):
    with open(f'{path}/model_config.json', 'r') as infile:
        config = json.load(infile, object_hook=get_decoder_hook())
    return config
