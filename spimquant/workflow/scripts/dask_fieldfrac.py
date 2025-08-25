from zarrnii import ZarrNii

# coiled not needed for this, since usually creating fieldfrac at a downsampled level
# if snakemake.config['use_coiled']:
#    from coiled import Cluster
#    cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=[4,30],idle_timeout='1 hour')
#    client = cluster.get_client()

in_orient = snakemake.config["in_orientation"]
orient_opt = {} if in_orient == None else {"orientation": in_orient}

ds_level = int(snakemake.wildcards.dslevel) - int(
    snakemake.config["segmentation_level"]
)

# this will downsample automatically based on the level
znimg_density_ds = ZarrNii.from_ome_zarr(
    snakemake.params.mask_uri, **orient_opt, level=ds_level
)

print(znimg_density_ds.darr)
znimg_density_ds.to_nifti(snakemake.output.fieldfrac_nii)
