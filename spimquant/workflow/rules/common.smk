from functools import cached_property
from snakebids import generate_inputs


class CachedConfig:
    """
    For lazy initialization of some variables.
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
