import nibabel as nib
from  zarrnii import ZarrNii, Transform
from dask.diagnostics import ProgressBar

#get channel index from omero metadata
zi = zarr.open(snakemake.input.ome_zarr)
attrs = zi['/'].attrs.asdict()
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)


#member function of floating image
flo_znimg = ZarrNii.from_path(snakemake.input.ome_zarr, channels=[channel_index])
ref_znimg = ZarrNii.from_path_as_ref(snakemake.input.ref_nii, channels=[channel_index],**snakemake.params.ref_opts)

out_znimg = flo_znimg.apply_transform(Transform.affine_ras_from_txt(snakemake.input.xfm_ras),ref_znimg=ref_znimg)

with ProgressBar():

    out_znimg.to_nifti(snakemake.output.nii,scheduler='single-threaded')



