from functools import cached_property
from snakebids import generate_inputs


class CachedConfig:
    """
    For lazy initialization of some variables. For example, the metadata read for a BIDS config needs to be done once
    and once only for rules that requires to read the BIDS dataset; and it does not need to be read at all for rules
    that don't require reading the BIDS dataset. This motivates the following pattern:

    class C:
        ...
        def get_inputs(self):
            if self.input is None:
                inputs = read_bids(...)
            return inputs
        ...

    which is more conveniently done using the @cached_property decorator as suggested by
    https://stackoverflow.com/questions/15226721/python-class-member-lazy-initialization
    So this class uses this method to do lazy intialization of variables
    """
    def __init__(self, config):
        self.config = config

    @cached_property
    def inputs(self):
        # Get input wildcards
        inputs = generate_inputs(
            bids_dir=self.config["bids_dir"],
            pybids_inputs=self.config["pybids_inputs"],
            pybidsdb_dir=self.config.get("pybidsdb_dir"),
            pybidsdb_reset=self.config.get("pybidsdb_reset"),
            derivatives=self.config.get("derivatives",None),
            participant_label=self.config.get("participant_label",None),
            exclude_participant_label=self.config.get("exclude_participant_label",None),
            validate=not self.config.get("plugins.validator.skip",False),
        )
        return inputs



def bids_tpl(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
    return str(Path(bids(root=root, tpl=template)) / bids(tpl=template, **entities))
