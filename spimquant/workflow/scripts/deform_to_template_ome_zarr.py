import nibabel as nib
from zarrnii import ZarrNii, AffineTransform, DisplacementTransform
from dask.diagnostics import ProgressBar

flo_znimg = ZarrNii.from_ome_zarr(
    snakemake.input.ome_zarr, channels=[snakemake.params.channel_index]
)
ref_znimg = ZarrNii.from_ome_zarr(
    snakemake.input.ref_nii,
    channels=[snakemake.params.channel_index],
    **snakemake.params.ref_opts,
    as_ref=True,
)


out_znimg = flo_znimg.apply_transform(
    DisplacementTransform.from_nifti(snakemake.input.warp_nii),
    AffineTransform.from_txt(snakemake.input.xfm_ras),
    ref_znimg=ref_znimg,
)

out_znimg.to_ome_zarr(snakemake.output.ome_zarr)
