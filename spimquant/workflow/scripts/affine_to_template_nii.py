import nibabel as nib
from  zarrnii import ZarrNii, AffineTransform
from dask.diagnostics import ProgressBar

#get channel index from omero metadata
zi = zarr.open(snakemake.input.ome_zarr)
attrs = zi['/'].attrs.asdict()
channel_labels = [channel_dict['label'] for channel_dict in attrs['omero']['channels']]
channel_index = channel_labels.index(snakemake.wildcards.stain)

in_orient = snakemake.config['in_orientation']
orient_opt = {} if in_orient == None else {'orientation': in_orient}


#member function of floating image
flo_znimg = ZarrNii.from_ome_zarr(snakemake.input.ome_zarr, channels=[channel_index],**orient_opt)
ref_znimg = ZarrNii.from_nifti(snakemake.input.ref_nii, channels=[channel_index],**snakemake.params.ref_opts,as_ref=True)

out_znimg = flo_znimg.apply_transform(AffineTransform.from_txt(snakemake.input.xfm_ras),ref_znimg=ref_znimg)

with ProgressBar():

    out_znimg.to_nifti(snakemake.output.nii,scheduler='single-threaded')



