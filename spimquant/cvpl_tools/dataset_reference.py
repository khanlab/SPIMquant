from datetime import datetime
from spimquant.cvpl_tools.array_key_dict import ArrayKeyDict
import spimquant.cvpl_tools.fs as fs
import numpy as np
import dask.array as da
import json


class DatapointReference:
    """
    Represents a pointer to location of data.
    The data is either an image or a chunk (data point) in the dataset.
    """

    def __init__(self, data_ref: str | list):
        assert isinstance(data_ref, (str, list)), 'ERROR: Datapoint reference should be either a string or a list!'
        self.data_ref = data_ref

    def to_json(self):
        return json.dumps(self.data_ref)

    def __str__(self):
        return self.to_json()

    def ref(self):
        return self.data_ref

    def has_multiple_images(self):
        return isinstance(self.data_ref, list)

    @classmethod
    def from_json(cls, string):
        return DatapointReference(json.loads(string))

    def read_as_np(self, read_setting: fs.ImReadSetting) -> np.ndarray:
        return fs.ImIO.read_single_image(read_setting, self.data_ref)


class DatasetReference:
    """
    This class' instance represents a dataset, and whose attributes include information on what data is in the dataset
    and when and how this reference is created.
    """

    def __init__(self, datapoint_refs: ArrayKeyDict[str, DatapointReference],
                 dataset_name,
                 creation_info):
        # create an empty dataset
        self.datapoint_refs = datapoint_refs
        self.creation_date = datetime.now()
        self.modification_info = creation_info
        self.name = dataset_name

    @classmethod
    def empty(cls):
        return DatasetReference(ArrayKeyDict(),
                                'Empty Dataset',
                                'This dataset is created as an empty dataset.')

    @classmethod
    def from_json(cls, json_str):
        d = json.loads(json_str)
        pairs = d['dataset_refs']
        refs = ArrayKeyDict((key, DatapointReference.from_json(datapt_ref_json_str))
                            for key, datapt_ref_json_str in pairs)
        dataset_ref = DatasetReference(refs, d['name'], d['modification_info'])
        t = d['creation_date']
        t = map(int, t[t.find('(') + 1: t.rfind(')')].split(','))
        dataset_ref.creation_date = datetime(*t)
        return dataset_ref

    @classmethod
    def load_from_json(cls, path):
        with open(path, 'r') as infile:
            json_str = infile.read()
        return DatasetReference.from_json(json_str)

    def to_json(self):
        refs = self.datapoint_refs
        pairs = {key: refs[key].to_json() for key in refs}
        t = repr(self.creation_date)
        d = {
            'name': self.name,
            'modification_info': self.modification_info,
            'creation_date': t,
            'dataset_refs': pairs,
        }
        return json.dumps(d)

    def save_as_json(self, path):
        with open(path, 'w') as outfile:
            outfile.write(self.to_json())


if __name__ == '__main__':
    in1 = 'C:/path/to/file1'
    r1 = DatapointReference(in1)
    in2 = ['C:/path/to/file2', 'C:/path/to/file3']
    r2 = DatapointReference(in2)
    assert str(r1) == '"C:/path/to/file1"'  # note the addition of double quotes
    assert DatapointReference.from_json(str(r1)).ref() == in1
    assert str(r2) == '["C:/path/to/file2", "C:/path/to/file3"]'
    assert DatapointReference.from_json(str(r2)).ref() == in2
