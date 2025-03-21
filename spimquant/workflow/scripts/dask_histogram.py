import numpy as np
from zarrnii import ZarrNii
from lib.utils import get_zarr_store
from dask.array import histogram

if snakemake.config['use_coiled']:
    from coiled import Cluster
    cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=30,idle_timeout='1 hour')
    client = cluster.get_client()

store = get_zarr_store(snakemake.params.spim_n4_uri)

in_orient = snakemake.config['in_orientation']
orient_opt = {} if in_orient == None else {'orientation': in_orient}


#we use the default level=0, since we are reading in the n4 output, which is already downsampled if level was >0
znimg_hires = ZarrNii.from_ome_zarr(store)

#note: histogram can take a weights dask array, which would be the perfect application for an upsampled brain mask

hist, bin_edges = histogram(znimg_hires.darr,**snakemake.params.histogram_opts)

hist.to_zarr(snakemake.params.histogram_uri)


