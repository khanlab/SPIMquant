from datetime import datetime
from .array_key_dict import ArrayKeyDict
import spimquant.cvpl_tools.fs as fs
import numpy as np
import dask.array as da
import json
from dataclasses import dataclass
import spimquant.cvpl_tools.strenc as strenc


@dataclass
class DatapointReference:
    """
    Represents a pointer to location of data.
    The data is either an image or a chunk (data point) in the dataset.
    """
    data_ref: str | list

    def __init__(self, data_ref: str | list):
        assert isinstance(data_ref, (str, list)), 'ERROR: Datapoint reference should be either a string or a list!'
        self.data_ref = data_ref

    def __str__(self):
        return str(self.data_ref)

    def ref(self):
        return self.data_ref

    def has_multiple_images(self):
        return isinstance(self.data_ref, list)

    def read_as_np(self, read_setting: fs.ImReadSetting) -> np.ndarray:
        return fs.ImIO.read_single_image(read_setting, self.data_ref)


@dataclass
class DatasetReference:
    """
    This class' instance represents a dataset, and whose attributes include information on what data is in the dataset
    and when and how this reference is created.
    """
    datapoint_refs: ArrayKeyDict[str, DatapointReference]
    creation_date: datetime
    creation_info: str  # a description of how this dataset is created
    name: str  # name of the dataset reference

    def __init__(self, datapoint_refs: ArrayKeyDict[str, DatapointReference],
                 dataset_name,
                 creation_info):
        # create an empty dataset
        self.datapoint_refs = datapoint_refs
        self.creation_date = datetime.now()
        self.creation_info = creation_info
        self.name = dataset_name

    @classmethod
    def empty(cls):
        return DatasetReference(ArrayKeyDict(),
                                'Empty Dataset',
                                'This dataset is created as an empty dataset.')


def test():
    in1 = 'C:/path/to/file1'
    r1 = DatapointReference(in1)
    in2 = ['C:/path/to/file2', 'C:/path/to/file3']
    r2 = DatapointReference(in2)
    assert str(r1) == 'C:/path/to/file1'
    assert len(r2.ref()) == 2 and r2.ref()[0] == 'C:/path/to/file2'
