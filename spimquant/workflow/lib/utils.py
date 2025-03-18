import zarr
import os
from upath import UPath as Path

def is_remote(uri_string):
    uri = Path(uri_string)
    if uri.protocol == 'gcs' or uri.protocol == 's3':
        return True
    else:
        return False

def is_remote_gcs(uri_string):
    uri = Path(uri_string)
    if uri.protocol == 'gcs':
        return True
    else:
        return False


def get_fsspec(uri_string,storage_provider_settings=None,creds=None):
    uri = Path(uri_string)
    if uri.protocol == 'gcs':
        from gcsfs import GCSFileSystem
        gcsfs_opts={}
        #gcsfs_opts={'project': storage_provider_settings['gcs'].get_settings().project,
        gcsfs_opts={'token': creds}
        fs = GCSFileSystem(**gcsfs_opts)
    elif uri.protocol == 's3':
        from s3fs import S3FileSystem
        s3fs_opts={'anon': False}
        fs = S3FileSystem(**s3fs_opts)
    elif uri.protocol == 'file' or uri.protocol == 'local' or uri.protocol == '':
        #assumed to be local file
        from fsspec.implementations.local import LocalFileSystem
        fs = LocalFileSystem()
    else:
        print(f'unsupported protocol for remote data')
    
    return fs


def get_storage_creds(uri,remote_creds):
    """for rules that deal with remote storage directly"""
    protocol = Path(uri).protocol
    if protocol == "gcs":
        # currently only works with gcs
        creds = os.path.expanduser(remote_creds) 
        return {"creds": creds}
    else:
        return {}




def get_channel_names(store):

    group = zarr.open_group(store, mode='r')

    # Access the .zattrs metadata as a dictionary
    metadata = dict(group.attrs)

    # Retrieve channel names
    return [chan['label'] for chan in metadata['multiscales'][0]['metadata']['omero']['channels']]

def get_channel_index(store, stain):
    channel_names = get_channel_names(store)
    return channel_names.index(stain)
    
def get_zarr_store(uri,creds=None):

    if is_remote(uri):
        #fs_args={'storage_provider_settings':snakemake.params.storage_provider_settings,'creds':snakemake.input.creds}
        fs_args={'creds':creds}
    else:
        fs_args={}

    fs = get_fsspec(uri,**fs_args)

    if Path(uri).suffix == '.zip':
        store = zarr.storage.ZipStore(Path(uri).path,dimension_separator='/',mode='r')
    else:
        store = zarr.storage.FSStore(Path(uri).path,fs=fs,dimension_separator='/',mode='r')

    return store

