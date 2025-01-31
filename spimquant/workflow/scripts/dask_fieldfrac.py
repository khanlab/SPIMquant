from zarrnii import ZarrNii

if snakemake.config['use_coiled']:
    from coiled import Cluster
    cluster = Cluster(name='coiled-snakemake',package_sync_ignore=['spimquant'],n_workers=[4,30],idle_timeout='1 hour')
    client = cluster.get_client()

in_orient = snakemake.config['in_orientation']
orient_opt = {} if in_orient == None else {'orientation': in_orient}

ds_level=int(snakemake.wildcards.dslevel)-int(snakemake.config["segmentation_level"])

znimg_mask = ZarrNii.from_ome_zarr(snakemake.params.mask_uri,**orient_opt)
znimg_density_ds = znimg_mask.downsample(level=ds_level,z_level_offset=0)
print(znimg_density_ds.darr)
znimg_density_ds.to_nifti(snakemake.output.fieldfrac_nii)
