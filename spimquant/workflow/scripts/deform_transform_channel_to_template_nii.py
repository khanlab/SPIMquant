import nibabel as nib
from  ome_zarr_neuro.transform import apply_transform, TransformSpec
from dask.diagnostics import ProgressBar

darr_interp = apply_transform(flo_img_path=snakemake.input.ome_zarr,
                ref_img_path=snakemake.input.ref_nii,
                transform_specs=[TransformSpec.displacement_from_nifti(snakemake.input.warp_nii),
                    TransformSpec.affine_ras_from_txt(snakemake.input.xfm_ras)],
                flo_channel=snakemake.params.channel_index,
                ref_chunks=(100,100,100))

ref_nib = nib.load(snakemake.input.ref_nii)


with ProgressBar():
    interp_vol = darr_interp.compute(scheduler='single-threaded') # uses too much memory otherwise
    #interp_vol = darr_interp.compute() # uses too much memory otherwise

out_nib = nib.Nifti1Image(interp_vol,
                       affine=ref_nib.affine,
                        header=ref_nib.header)


out_nib.to_filename(snakemake.output.nii)



