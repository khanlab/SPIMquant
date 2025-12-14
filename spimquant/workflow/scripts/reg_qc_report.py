#!/usr/bin/env python
# coding: utf-8

# In[1]:


######## snakemake preamble start (automatically inserted, do not edit) ########
######## snakemake preamble end #########


# # Registration Quality Control Report
# 
# This notebook provides comprehensive quality control visualizations for template registration,
# including comparisons of template and transformed subject images, overlay visualizations,
# label maps, and warp field visualizations.

# In[2]:


import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from nilearn import plotting
from pathlib import Path


# In[3]:


# Load inputs from snakemake object
template_path = snakemake.input.template
subject_path = snakemake.input.subject
warped_affine_path = snakemake.input.warped_affine
warped_deform_path = snakemake.input.warped_deform
warp_field_path = snakemake.input.warp
dseg_path = snakemake.input.dseg

output_notebook = snakemake.output.notebook

# Get wildcards for report title
subject = snakemake.wildcards.subject
template = snakemake.wildcards.template
stain = snakemake.wildcards.stain

print(f"Processing registration QC for subject {subject}, stain {stain}, template {template}")


# In[4]:


# Load NIfTI images
template_img = nib.load(template_path)
subject_img = nib.load(subject_path)
warped_affine_img = nib.load(warped_affine_path)
warped_deform_img = nib.load(warped_deform_path)
warp_img = nib.load(warp_field_path)
dseg_img = nib.load(dseg_path)

print(f"Template shape: {template_img.shape}")
print(f"Subject shape: {subject_img.shape}")
print(f"Warped (affine) shape: {warped_affine_img.shape}")
print(f"Warped (deform) shape: {warped_deform_img.shape}")
print(f"Warp field shape: {warp_img.shape}")
print(f"Segmentation shape: {dseg_img.shape}")


# ## 1. Affine Registration Quality
# 
# Comparison of template and affine-transformed subject

# In[44]:


# Create figure for affine registration comparison
fig, axes = plt.subplots(2, 1, figsize=(15, 10))

fig.suptitle(f"Affine Registration: {subject} to {template}", fontsize=16)

# Template orthoviews
display = plotting.plot_anat(template_img, title="Template", 
                              display_mode='ortho', figure=fig, axes=axes[0],cut_coords=(0,0,0))

# Affine warped subject orthoviews
display = plotting.plot_anat(warped_affine_img, title="Affine Warped Subject",
                              display_mode='ortho', figure=fig, axes=axes[1],cut_coords=(0,0,0),dim=-1.5)

plt.tight_layout()
plt.show()


# In[18]:


# Overlay visualization for affine registration
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Affine Registration Overlay: Template (red) + Subject (blue/green)", fontsize=16)

display = plotting.plot_anat(template_img, title="Affine Overlay",
                              display_mode='z', cut_coords=7, figure=fig)
display.add_overlay(warped_affine_img, cmap='viridis', transparency=0.5)
plt.show()


# ## 2. Deformable Registration Quality
# 
# Comparison of template and deformably-transformed subject

# In[20]:


# Create figure for deformable registration comparison
fig, axes = plt.subplots(2, 1, figsize=(15, 10))
fig.suptitle(f"Deformable Registration: {subject} to {template}", fontsize=16)

# Template orthoviews
display = plotting.plot_anat(template_img, title="Template", 
                              display_mode='ortho', figure=fig, axes=axes[0])

# Deformably warped subject orthoviews
display = plotting.plot_anat(warped_deform_img, title="Deformably Warped Subject",
                              display_mode='ortho', figure=fig, axes=axes[1])

plt.tight_layout()
plt.show()


# In[21]:


# Overlay visualization for deformable registration
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Deformable Registration Overlay: Template (red) + Subject (blue/green)", fontsize=16)

display = plotting.plot_anat(template_img, title="Deformable Overlay",
                              display_mode='z', cut_coords=7, figure=fig)
display.add_overlay(warped_deform_img, cmap='viridis', transparency=0.5)
plt.show()


# ## 3. Edge/Contour Overlays
# 
# Visualizing registration quality with edge contours

# In[31]:


# Affine registration with contours
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Affine Registration with Contours", fontsize=16)

display = plotting.plot_anat(warped_affine_img, title="Affine with Template Contours",
                              display_mode='z', cut_coords=7, figure=fig,dim=-1)
display.add_contours(template_img, levels=[100, 1000,2000], colors='r', transparency=0.7)
plt.show()


# In[32]:


# Deformable registration with contours
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Deformable Registration with Contours", fontsize=16)

display = plotting.plot_anat(warped_deform_img, title="Deformable with Subject Contours",
                              display_mode='z', cut_coords=7, figure=fig,dim=-1)
display.add_contours(template_img, levels=[100, 1000,2000], colors='r', transparency=0.7)
plt.show()


# ## 4. Label Map Visualization
# 
# Visualization of atlas segmentation on template and transformed subject

# In[34]:


# Label map overlay on template
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Atlas Segmentation on Template", fontsize=16)

display = plotting.plot_roi(dseg_img, bg_img=template_img,
                             title="Atlas Labels on Template",
                             display_mode='z', cut_coords=7, 
                             cmap='tab20', figure=fig)
plt.show()


# In[35]:


# Label map overlay on warped subject (deformable)
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Atlas Segmentation on Warped Subject", fontsize=16)

display = plotting.plot_roi(dseg_img, bg_img=warped_deform_img,
                             title="Atlas Labels on Deformably Warped Subject",
                             display_mode='z', cut_coords=7,
                             cmap='tab20', figure=fig)
plt.show()


# ## 5. Warp Field Visualization
# 
# Visualization of the deformation field magnitude and components

# In[39]:


# Calculate warp field magnitude
warp_data = warp_img.get_fdata()
if len(warp_data.shape) == 5 and warp_data.shape[-1] == 3:
    # Calculate magnitude of displacement
    warp_magnitude = np.sqrt(np.sum(warp_data**2, axis=-1))

    # Create a new NIfTI image for the magnitude
    warp_mag_img = nib.Nifti1Image(warp_magnitude, warp_img.affine, warp_img.header)

    # Visualize warp magnitude
    fig = plt.figure(figsize=(20, 5))
    fig.suptitle(f"Deformation Field Magnitude", fontsize=16)

    display = plotting.plot_stat_map(warp_mag_img, bg_img=template_img,
                                      title="Warp Magnitude",
                                      display_mode='z', cut_coords=7,
                                      cmap='hot', figure=fig)
    plt.show()

    print(f"Warp magnitude statistics:")
    print(f"  Mean: {np.mean(warp_magnitude):.3f}")
    print(f"  Std:  {np.std(warp_magnitude):.3f}")
    print(f"  Max:  {np.max(warp_magnitude):.3f}")
else:
    print(f"Warp field has unexpected shape: {warp_data.shape}")


# In[40]:


# Visualize individual warp components
if len(warp_data.shape) == 5 and warp_data.shape[-1] == 3:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle(f"Deformation Field Components (X, Y, Z)", fontsize=16)

    component_names = ['X', 'Y', 'Z']
    for i, (ax, name) in enumerate(zip(axes, component_names)):
        component_img = nib.Nifti1Image(warp_data[..., i], warp_img.affine, warp_img.header)
        display = plotting.plot_stat_map(component_img, bg_img=template_img,
                                          title=f"Warp {name} component",
                                          display_mode='z', cut_coords=5,
                                          cmap='coolwarm', symmetric_cbar=True,
                                          axes=ax)
    plt.tight_layout()
    plt.show()


# ## Summary
# 
# Registration QC report completed successfully.

# In[38]:


print("="*80)
print(f"Registration QC Report Summary")
print("="*80)
print(f"Subject: {subject}")
print(f"Template: {template}")
print(f"Stain: {stain}")
print(f"Output notebook: {output_notebook}")
print("="*80)
print("QC visualizations generated:")
print("  ✓ Affine registration comparison")
print("  ✓ Affine registration overlay")
print("  ✓ Deformable registration comparison")
print("  ✓ Deformable registration overlay")
print("  ✓ Affine registration with contours")
print("  ✓ Deformable registration with contours")
print("  ✓ Atlas segmentation on template")
print("  ✓ Atlas segmentation on warped subject")
print("  ✓ Deformation field magnitude")
print("  ✓ Deformation field components")
print("="*80)


# In[ ]:




