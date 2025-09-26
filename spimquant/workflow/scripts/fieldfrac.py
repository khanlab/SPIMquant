from zarrnii import ZarrNii

ds_level = int(snakemake.wildcards.dslevel) - int(
    snakemake.config["segmentation_level"]
)

# this will downsample automatically based on the level
znimg_density_ds = ZarrNii.from_ome_zarr(
    snakemake.input.mask,
    level=ds_level,
    downsample_near_isotropic=True,
)

print(znimg_density_ds.darr)
znimg_density_ds.to_nifti(snakemake.output.fieldfrac_nii)
