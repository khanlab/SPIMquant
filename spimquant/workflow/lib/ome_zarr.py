import zarr

def get_channel_names(path):

    group = zarr.open_group(path, mode='r')

    # Access the .zattrs metadata as a dictionary
    metadata = dict(group.attrs)

    # Retrieve channel names
    return [chan['label'] for chan in metadata['multiscales'][0]['metadata']['omero']['channels']]

    return channel_names


