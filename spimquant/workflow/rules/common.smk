from upath import UPath
from pathlib import Path

def bids_tpl(root, template, **entities):
    """bids() wrapper for files in tpl-template folder"""
    return str(Path(bids(root=root, tpl=template)) / bids(tpl=template, **entities))

def get_storage_creds(uri):
    """for rules that deal with remote storage directly"""
    protocol = UPath(uri).protocol
    if protocol == "gcs":
        # currently only works with gcs
        creds = os.path.expanduser(config["remote_creds"])
        return {"creds": creds}
    else:
        return {}
