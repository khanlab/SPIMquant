from zarrnii import ZarrNii
from dask.diagnostics import ProgressBar
from lib.utils import get_zarr_store, get_channel_index


store = get_zarr_store(snakemake.params.in_zarr)

channel_index = get_channel_index(store, snakemake.wildcards.stain)

level=int(snakemake.wildcards.level)

in_orient = snakemake.config['in_orientation'] #TODO: this is a bit ugly - update the ZarrNii to recognize None  as an orientation, so we can just pass in_orientation
orient_opt = {} if in_orient == None else {'orientation': in_orient}

with ProgressBar():
    # ZarrNii now implicitly downsamples xy and z if warranted by chosen level
    ZarrNii.from_ome_zarr(store,level=level,channels=[channel_index],**orient_opt).to_nifti(snakemake.output.nii)

    
