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




