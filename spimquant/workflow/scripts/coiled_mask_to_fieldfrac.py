
# use the binary mask to produce a density map:


def get_downsampled_density(znimg,along_x=1,along_y=1,along_z=1):
    """ downsamples by local mean"""

    if znimg.axes_nifti:
        axes ={0:1,
                                                    1: along_x,
                                                    2: along_y,
                                                    3: along_z}
    else:
        axes ={0:1,
                                                    1: along_z,
                                                    2: along_y,
                                                    3: along_x}

    def norm_sum(x, *args, **kwargs):
        
        return np.sum(x, *args, **kwargs) / float(along_x*along_y*along_z)
        
    #coarsen performs a reduction in a local neighbourhood, defined by axes
    darr_scaled = da.coarsen(norm_sum,x=znimg.darr.astype('float32'),axes=axes, trim_excess=True) 


    #we need to also update the affine, scaling by the ds_factor
    scaling_matrix = np.diag((along_x,along_y,along_z,1))
    new_vox2ras = scaling_matrix @ znimg.vox2ras.affine

    return ZarrNii.from_darr(darr_scaled,vox2ras=new_vox2ras,axes_nifti=znimg.axes_nifti)



#this saves a downsampled nifti from the mask
#znimg_mask = ZarrNii.from_path(snakemake.params.mask_uri)
#znimg_density_ds = get_downsampled_density(znimg_mask,along_x=32,along_y=32,along_z=16)
#znimg_density_ds.to_nifti('test_density_level5iso_fieldfrac.nii')



