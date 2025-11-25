from zarrnii import ZarrNii

target_level = int(snakemake.wildcards.level)
hires_level = int(snakemake.params.hires_level)

downsampling_level = target_level - hires_level
if downsampling_level < 0:
    raise ValueError("Target level for fieldfrac is smaller than the input level!")


# this will downsample automatically based on the level
znimg_density_ds = ZarrNii.from_ome_zarr(
    snakemake.input.mask,
    level=downsampling_level,
    downsample_near_isotropic=True,
    **snakemake.params.zarrnii_kwargs,
)

print(znimg_density_ds.darr)
znimg_density_ds.to_nifti(snakemake.output.fieldfrac_nii)
