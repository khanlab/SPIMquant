import zarr
import nibabel as nib
from  zarrnii import ZarrNii, Transform
from dask.distributed import Client

#compare target resolution to input pyramid to determine what resolution to get 
zooms_mm = snakemake.params.ref_opts['zooms'][-1] #x zoom
print(f'zooms_mm: {zooms_mm}')

from lib.ome_zarr import get_metadata

metadata = get_metadata(snakemake.params.ome_zarr)
coordinate_transforms = metadata['coordinateTransformations']
print(coordinate_transforms)
level=0
for i,c in enumerate(coordinate_transforms):
    x_zoom = c[0]['scale'][-1] 
    if x_zoom < zooms_mm:
        print(f'x_zoom: {x_zoom} is smaller, so setting level={i}')
        level=i
    else:
        print(f'x_zoom: {x_zoom} is larger, so using prev level={level}')
        break

print(f'level: {level}')


client = Client(n_workers=6, threads_per_worker=2,processes=False)
print(client.dashboard_link)

#get channel index from omero metadata
zi = zarr.open(snakemake.params.ome_zarr)
attrs = zi['/'].attrs.asdict()
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)


#member function of floting image
flo_znimg = ZarrNii.from_path(snakemake.params.ome_zarr, channels=[channel_index], level=level,**snakemake.params.flo_opts)
ref_znimg = ZarrNii.from_path_as_ref(snakemake.input.ref_nii, channels=[channel_index],**snakemake.params.ref_opts)


if snakemake.params.do_downsample:
    flo_znimg = flo_znimg.downsample(**snakemake.params.downsample_opts)


deform_znimg = flo_znimg.apply_transform(Transform.displacement_from_nifti(snakemake.input.warp_nii),
                    Transform.affine_ras_from_txt(snakemake.input.xfm_ras),
                    ref_znimg=ref_znimg)


deform_znimg.to_nifti(snakemake.output.nii)



client.close()
