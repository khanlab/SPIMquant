import nibabel as nib
from  ome_zarr_neuro.transform import DaskImage, TransformSpec
from dask.diagnostics import ProgressBar


#member function of floting image
flo_dimg = DaskImage.from_path(snakemake.input.ome_zarr, channels=[snakemake.params.channel_index])
ref_dimg = DaskImage.from_path_as_ref(snakemake.input.ref_nii, channels=[snakemake.params.channel_index],chunks=snakemake.params.chunks,zooms=snakemake.params.zooms)

out_dimg = flo_dimg.apply_transform(TransformSpec.affine_ras_from_txt(snakemake.input.xfm_ras),ref_dimg=ref_dimg)

with ProgressBar():

    out_dimg.to_nifti(snakemake.output.nii,scheduler='single-threaded')



