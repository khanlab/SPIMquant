from coiled import Cluster
from zarrnii import ZarrNii
import dask.array as da
import numpy as np

cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=[4,20])
client = cluster.get_client()


ds_level=int(snakemake.wildcards.dslevel)

znimg_mask = ZarrNii.from_path(snakemake.params.mask_uri)
znimg_density_ds = znimg_mask.downsample(level=ds_level,do_normalized_sum=True)
print(znimg_density_ds.darr)
znimg_density_ds.to_nifti(snakemake.output.fieldfrac_nii)
