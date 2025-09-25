from lib.utils import get_zarr_store
from zarrnii import ZarrNii, BiasFieldCorrection
from dask.diagnostics import ProgressBar

store = get_zarr_store(snakemake.params.spim_uri)

hires_level = int(snakemake.wildcards.level)
ds_level = int(snakemake.wildcards.dslevel)

znimg = ZarrNii.from_ome_zarr(
    store,
    channel_labels=[snakemake.wildcards.stain],
    level=hires_level,
    downsample_near_isotropic=True,
)


print('compute bias field correction')
with ProgressBar():
    
    # Apply bias field correction
    znimg_corrected = znimg.apply_scaled_processing(
        BiasFieldCorrection(sigma=5.0),
        downsample_factor=2,
        use_temp_zarr=False #bug when set to True - look into zarrnii

    )

    # write to ome_zarr
    znimg_corrected.to_ome_zarr(snakemake.params.spim_n4_uri, max_layer=5)


    #write to nifti as test
    znimg.to_nifti('orig.nii')
    znimg_corrected.to_nifti('orig_n4.nii')
