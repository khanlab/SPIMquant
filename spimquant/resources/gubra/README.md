Gubra AutoF template 
Modified from: https://github.com/Gubra-ApS/LSFM-mouse-brain-atlas

(posted files did not have correct header information (voxel size), was not in similar physical space to ABAv3, and had a different label mapping. This corrects by using the below c3d commands for spacing reslicing, and a python script to remap the labels. Note: the rigid transform was found using itksnap (align by center, then MI rigid reg)

e.g.:
```
c3d gubra_template_olf.nii.gz -spacing 0.025x0.025x0.025mm -o gubra_template_olf_spacing.nii.gz
c3d ../../ABAv3/P56_Atlas.nii.gz  gubra_template_olf_spacing.nii.gz  -reslice-itk ../../gubra/gubra_withspacing_to_aba_rigid_itk.txt -o ../../gubra/gubra_template_olf_spacing_reslice.nii.gz
c3d -int 0 gubra_template_olf_spacing_reslice.nii.gz gubra_ano_olf_spacing_remap.nii.gz -reslice-itk gubra_withspacing_to_aba_rigid_itk.txt -o gubra_ano_olf_spacing_remap_reslice.nii.gz
```
