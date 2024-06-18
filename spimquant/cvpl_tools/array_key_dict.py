import os
import json


class ArrayKeyDict(dict):
    """
    This class ensures a dictionary with a customizable ordering imposed over the keys.
    Intended use cases additional to a Python dictionary:
    - Matching up two or more dictionaries' items by keys
    - Dictionary subsetting and merging
    - Modifying the key of an item while maintaining the position of the item in the dict
    - Looping over the keys in specific orders (may be different from the order they are inserted)

    Example:
        d = ArrayKeyDict([(3, 'a'), (4, 'b')])
        print(d.key_at(0))  # outputs 3
        for k in d:
            print(d[k], end='')  # outputs 'ab'
        print()
        d2 = ArrayKeyDict([(3, 'c'), (5, 'e')])
        comm_keys = d.common_keylist_with(d2)  # [3]
        d3 = d.subset_dict(comm_keys)  # now becomes [(3, 'a')]
        for k in d3:
            print(d3[k], end='')  # outputs 'a'
    """
    def __init__(self, pairs=()):
        """
        Create a dictionary; if the pairs parameter is provided, create a dictionary with items following the
        orders in the provided pairs
        :param pairs: An iterable of key, value pairs that specifies the initial items in the dictionary
        """
        super().__init__()
        self._key_arr = []

        # ensure ordering
        for k, v in pairs:
            self[k] = v

    def __setitem__(self, key, value):
        """
        Ensure the newly inserted value is at the end of key array;
        If the key already exists, do not change its order in the array
        """
        flag = key in self
        super().__setitem__(key, value)
        if not flag:
            self._key_arr.append(key)

    def keys(self) -> list:
        return list(self._key_arr)

    def values(self) -> list:
        return [self[k] for k in self._key_arr]

    def key_at(self, idx):
        return self._key_arr[idx]

    def __iter__(self):
        return iter(self._key_arr)

    def to_json(self):
        pairs = {key: self[key] for key in self}
        return json.dumps(pairs)

    @classmethod
    def from_json(cls, json_str):
        pairs = json.loads(json_str)
        return ArrayKeyDict((key, pairs[key]) for key in pairs)

    def pop(self, key, default=None):
        key_in_arr = key in self._key_arr
        if not key_in_arr:
            if default is None:
                raise KeyError(f'ERROR: Key {key} not found in array when pop is called!')
            else:
                return default

        self._key_arr.remove(key)
        item = super().pop(key)
        return item

    def __delitem__(self, key):
        key_in_arr = key in self._key_arr
        if not key_in_arr:
            raise KeyError(f'ERROR: Key {key} not found in array when del is called!')

        self._key_arr.remove(key)
        super().__delitem__(key)

    def rename_key(self, old_key, new_key):
        if new_key == old_key:
            return
        assert new_key not in self, f'ERROR: Attempt to rename to an existing key {new_key}'
        super().__setitem__(new_key, super().pop(old_key))
        for i in range(len(self._key_arr)):
            if self._key_arr[i] == old_key:
                self._key_arr[i] = new_key
                break

    def reorder_keys(self, ordered_key_arr):
        # no new keys & length equal
        assert len(self._key_arr) == len(ordered_key_arr), f'ERROR: lengths of key arrays are different during reorder!'\
                                                           f' old {len(self._key_arr)} and new {len(ordered_key_arr)}'
        # no repetition
        assert len(ordered_key_arr) == len(set(ordered_key_arr)), f'ERROR: Repetition found in input {ordered_key_arr}!'
        # all keys are in the original dict
        for key in ordered_key_arr:
            assert key in self, f'ERROR: New key array contains value not present in the original key array: {key}'

        # all checking passed, we can safely replace the array
        self._key_arr = ordered_key_arr

    def subset_dict(self, subset_key_arr):
        # no new keys & length smaller
        assert len(self._key_arr) >= len(subset_key_arr), f'ERROR: subset key arr has larger length {len(subset_key_arr)}'\
                                                            f' than the original {len(self._key_arr)}'
        # no repetition
        assert len(subset_key_arr) == len(set(subset_key_arr)), f'ERROR: Repetition found in input {subset_key_arr}!'
        # all keys are in the original dict
        for key in subset_key_arr:
            assert key in self, f'ERROR: New key array contains value not present in the original key array: {key}'
        pairs = ((key, self[key]) for key in subset_key_arr)
        return ArrayKeyDict(pairs)

    def clear(self):
        self._key_arr.clear()
        super().clear()

    def trim_key_prefix_and_suffix(self, trim_front=0, trim_end=0):
        for i in range(len(self)):
            old_key = self._key_arr[i]
            new_key = old_key[trim_front:len(old_key) - trim_end]
            self.rename_key(old_key, new_key)

    def trim_key_common_prefix_and_suffix(self):
        keys = [s for s in self._key_arr]
        rkeys = [s[::-1] for s in self._key_arr]
        cprefix = os.path.commonprefix(keys)
        csuffix = os.path.commonprefix(rkeys)
        self.trim_key_prefix_and_suffix(len(cprefix), len(csuffix))

    def common_keylist_with(self, dicts):
        """
        return a Python list containing the common keys of this dictionary (self) and the provided other dictionary.
        The order of which the keys will appear in the result is the same as they are in self.
        :param dicts: either an ArrayKeyDict or a list of ArrayKeyDict
        """
        if type(dicts) is ArrayKeyDict:
            dicts = [dicts]
        keylist = []
        for key in self:
            is_common = True
            for d in dicts:
                if key not in d:
                    is_common = False
                    break
            if is_common:
                keylist.append(key)
        return keylist


def test():
    test_dict1 = ArrayKeyDict()
    test_dict2 = ArrayKeyDict([['file1', 1], ['file2', 3], ['file3', 2]])
    assert len(test_dict1) == 0 and len(test_dict2) == 3
    test_dict2['file2'] = 4
    assert len(test_dict2) == 3 and test_dict2['file2'] == 4
    test_dict2['file4'] = 2
    assert test_dict2.keys() == ['file1', 'file2', 'file3', 'file4'] and len(test_dict2) == 4
    assert test_dict2.key_at(2) == 'file3'

    for key in reversed(test_dict2):
        test_dict1[key] = test_dict2[key]
    assert test_dict1.keys() == ['file4', 'file3', 'file2', 'file1']
    test_dict1.rename_key('file1', 'f1')
    assert test_dict1.key_at(3) == 'f1' and test_dict1['f1'] == 1
    test_dict1.reorder_keys(['f1', 'file4', 'file2', 'file3'])
    assert test_dict1.key_at(1) == 'file4' and test_dict1.keys() == ['f1', 'file4', 'file2', 'file3']

    # prevent regression where rename_key() incorrectly push the renamed key to end of array
    test_dict1.rename_key('file2', 'f2')
    assert len(test_dict1) == 4  # prevent regression when rename_key adds length to _key_arr
    assert test_dict1.key_at(2) == 'f2'

    test_dict3 = ArrayKeyDict([['file1', 1], ['file2', 3], ['file3', 2]])
    assert test_dict1.common_keylist_with([test_dict2, test_dict3]) == ['file3']
    assert test_dict3.subset_dict(['file3', 'file1']).keys() == ['file3', 'file1']

    assert test_dict1.keys() == ['f1', 'file4', 'f2', 'file3']
    test_dict1.trim_key_common_prefix_and_suffix()
    assert test_dict1.key_at(0) == '1' and test_dict1.key_at(1) == 'ile4'

    test_dict1.clear()
    assert len(test_dict1) == 0 and len(test_dict2) == 4
    test_dict2.clear()
    assert len(test_dict1) == 0 and len(test_dict2) == 0

    print('All tests successfully completed')