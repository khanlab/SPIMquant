import ome_zarr
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader


def get_metadata(path):

    # Parse the URL and load the OME-Zarr file
    store = parse_url(path)
    reader = Reader(store)
    nodes = list(reader())

    # Assume the first node is the image you want to read metadata from
    image_node = nodes[0]

    # Access the metadata
    metadata = image_node.metadata
    return metadata 

def get_channel_names(path):

    # Parse the URL and load the OME-Zarr file
    store = parse_url(path)
    reader = Reader(store)
    nodes = list(reader())

    # Assume the first node is the image you want to read metadata from
    image_node = nodes[0]

    # Access the metadata
    metadata = image_node.metadata

    # Retrieve channel names
    channel_names = metadata['channel_names']

    return channel_names


