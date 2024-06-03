import nibabel as nib
from  zarrnii import ZarrNii, Transform
from dask.diagnostics import ProgressBar

flo_znimg = ZarrNii.from_path(snakemake.input.ome_zarr, channels=[snakemake.params.channel_index])
ref_znimg = ZarrNii.from_path_as_ref(snakemake.input.ref_nii, channels=[snakemake.params.channel_index],**snakemake.params.ref_opts)


out_znimg = flo_znimg.apply_transform(Transform.displacement_from_nifti(snakemake.input.warp_nii),
                                Transform.affine_ras_from_txt(snakemake.input.xfm_ras),ref_znimg=ref_znimg)

out_znimg.to_ome_zarr(snakemake.output.ome_zarr)

