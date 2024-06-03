import zarr
import nibabel as nib
from  zarrnii import ZarrNii, Transform
from dask.distributed import Client

client = Client(n_workers=4, threads_per_worker=2,processes=False)
print(client.dashboard_link)

#get channel index from omero metadata
zi = zarr.open(snakemake.input.ome_zarr)
attrs = zi['/'].attrs.asdict()
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)


#member function of floting image
flo_znimg = ZarrNii.from_path(snakemake.input.ome_zarr, channels=[channel_index], **snakemake.params.flo_opts)
ref_znimg = ZarrNii.from_path_as_ref(snakemake.input.ref_nii, channels=[channel_index],**snakemake.params.ref_opts)

if snakemake.params.do_downsample:
    flo_ds_znimg = flo_znimg.downsample(**snakemake.params.downsample_opts)


deform_znimg = flo_ds_znimg.apply_transform(Transform.displacement_from_nifti(snakemake.input.warp_nii),
                    Transform.affine_ras_from_txt(snakemake.input.xfm_ras),
                    ref_znimg=ref_znimg)


deform_znimg.to_nifti(snakemake.output.nii)



client.close()
