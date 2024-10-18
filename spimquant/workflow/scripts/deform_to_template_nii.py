import zarr
import nibabel as nib
from  zarrnii import ZarrNii, Transform
import dask
from lib.cloud_io import get_fsspec, is_remote
from upath import UPath as Path
from dask.diagnostics import ProgressBar

in_zarr_url = snakemake.params.ome_zarr

dask.config.set(scheduler='threads', num_workers=snakemake.threads)  

if is_remote(in_zarr_url):
    #fs_args={'storage_provider_settings':snakemake.params.storage_provider_settings,'creds':snakemake.input.creds}
    fs_args={'creds':snakemake.input.creds}
else:
    fs_args={}

fs = get_fsspec(in_zarr_url,**fs_args)

if Path(in_zarr_url).suffix == '.zip':
    store = zarr.storage.ZipStore(Path(in_zarr_url).path,dimension_separator='/',mode='r')
else:
    store = zarr.storage.FSStore(Path(in_zarr_url).path,fs=fs,dimension_separator='/',mode='r')


zi = zarr.open(store=store,mode='r')


#get channel index from omero metadata
attrs = zi['/'].attrs.asdict()
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)


#member function of floting image
flo_znimg = ZarrNii.from_path(store, channels=[channel_index], **snakemake.params.flo_opts)
ref_znimg = ZarrNii.from_path_as_ref(snakemake.input.ref_nii, channels=[channel_index],**snakemake.params.ref_opts)

if snakemake.params.do_downsample:
    flo_ds_znimg = flo_znimg.downsample(**snakemake.params.downsample_opts)


deform_znimg = flo_ds_znimg.apply_transform(Transform.displacement_from_nifti(snakemake.input.warp_nii),
                    Transform.affine_ras_from_txt(snakemake.input.xfm_ras),
                    ref_znimg=ref_znimg)


with ProgressBar():
    deform_znimg.to_nifti(snakemake.output.nii)



