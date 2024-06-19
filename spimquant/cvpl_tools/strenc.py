"""
This file provides a JSON encoder and its corresponding decoder ('s hook).
To encode an object supported by this encoder:
json_string = json.dumps(object, cls=strenc.get_encoder())
To decode an object:
object = json.loads(json_string, object_hook=strenc.get_decoder_hook())

Note the behavior of json encode/decode will encode a subfield that is a dataclass as plain dict without
'__type__' field, as opposed to directly encoding instance of a dataclass class.
"""


import json
import dataclasses
from datetime import datetime


_Encoder = None
_decoder_hook = None


def get_encoder():
    if _Encoder is None:
        init()
    return _Encoder


def get_decoder_hook():
    if _decoder_hook is None:
        init()
    return _decoder_hook


def init():
    from .fs import ImReadSetting, ImWriteSetting
    from .array_key_dict import ArrayKeyDict
    from .dataset_reference import DatapointReference, DatasetReference


    # in order to support a more automatic conversion of dataclasses
    dataclasses_names = (
        ImReadSetting, ImWriteSetting, DatapointReference, DatasetReference
    )
    # here we use .__name__ to get class names without prefixes to ensure consistency of the class names
    # across Python invokations
    str_to_cls_map = {cls.__name__: cls for cls in dataclasses_names}


    class EncoderDef(json.JSONEncoder):
        """
        Turns an object to string and back; support dataclasses and datetime in addition to
        common Python types.
        """
        def default(self, o):
            if isinstance(o, ArrayKeyDict):
                d = {key: o[key] for key in o}
                d['__type__'] = 'array_key_dict.ArrayKeyDict'
                return d
            elif isinstance(o, datetime):
                d = {
                    '__type__': 'datetime.datetime',
                    'v1': o.isoformat()
                }
                return d
            elif dataclasses.is_dataclass(o):
                # reference: https://stackoverflow.com/questions/51286748/make-the-python-json-encoder-support-pythons-new-dataclasses
                d = dataclasses.asdict(o)
                d['__type__'] = type(o).__name__
                return d
            return super().default(o)


    def decoder_hook_def(o):
        if '__type__' in o:
            ty = o.pop('__type__')
            if ty == 'array_key_dict.ArrayKeyDict':
                return ArrayKeyDict((key, o[key]) for key in o)
            elif ty == 'datetime.datetime':
                return datetime.fromisoformat(o['v1'])
            else:
                cls = str_to_cls_map[ty]
                return cls(**o)
        return o


    global _Encoder, _decoder_hook
    _Encoder = EncoderDef
    _decoder_hook = decoder_hook_def


def test():
    from .fs import ImReadSetting
    from .dataset_reference import DatasetReference
    set = ImReadSetting()
    set_json = json.dumps(set, cls=get_encoder())
    set2: ImReadSetting = json.loads(set_json, object_hook=get_decoder_hook())
    assert set2.im_format == set.im_format
    dr = DatasetReference.empty()
    dr.im_read_setting = set
    dr_json = json.dumps(dr, cls=get_encoder())
    dr2: DatasetReference = json.loads(dr_json, object_hook=get_decoder_hook())
    assert dr.name == dr2.name
