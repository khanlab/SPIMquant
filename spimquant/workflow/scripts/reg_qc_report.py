#!/usr/bin/env python
# coding: utf-8

"""
Registration Quality Control Report Generator

This script generates an HTML report with comprehensive quality control 
visualizations for template registration, including comparisons of template 
and transformed subject images, overlay visualizations, label maps, and 
warp field visualizations.
"""

import base64
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from nilearn import plotting

# Load inputs from snakemake object
template_path = snakemake.input.template
subject_path = snakemake.input.subject
warped_affine_path = snakemake.input.warped_affine
warped_deform_path = snakemake.input.warped_deform
warp_field_path = snakemake.input.warp
dseg_path = snakemake.input.dseg

output_html = snakemake.output.report_html

# Get wildcards for report title
subject = snakemake.wildcards.subject
template = snakemake.wildcards.template
stain = snakemake.wildcards.stain

print(
    f"Processing registration QC for subject {subject}, stain {stain}, template {template}"
)

# HTML components storage
html_figures = []


def fig_to_base64(fig):
    """Convert matplotlib figure to base64 encoded string for HTML embedding."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return f"data:image/png;base64,{img_base64}"


def add_figure_to_html(fig, caption=""):
    """Add a matplotlib figure to the HTML report."""
    img_data = fig_to_base64(fig)
    html_figures.append({"image": img_data, "caption": caption})


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

# Store image info for report
image_info = {
    "template_shape": template_img.shape,
    "subject_shape": subject_img.shape,
    "warped_affine_shape": warped_affine_img.shape,
    "warped_deform_shape": warped_deform_img.shape,
    "warp_shape": warp_img.shape,
    "dseg_shape": dseg_img.shape,
}

# Section 1: Affine Registration Quality
print("Generating affine registration comparison...")
fig, axes = plt.subplots(2, 1, figsize=(15, 10))
fig.suptitle(f"Affine Registration: {subject} to {template}", fontsize=16)

display = plotting.plot_anat(
    template_img,
    title="Template",
    display_mode="ortho",
    figure=fig,
    axes=axes[0],
    cut_coords=(0, 0, 0),
)
display = plotting.plot_anat(
    warped_affine_img,
    title="Affine Warped Subject",
    display_mode="ortho",
    figure=fig,
    axes=axes[1],
    cut_coords=(0, 0, 0),
    dim=-1.5,
)

plt.tight_layout()
add_figure_to_html(
    fig,
    "Affine registration comparison: template (top) vs. affine-transformed subject (bottom)",
)

print("Generating affine overlay...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Affine Registration Overlay: Template + Subject", fontsize=16)

display = plotting.plot_anat(
    template_img, title="Affine Overlay", display_mode="z", cut_coords=7, figure=fig
)
display.add_overlay(warped_affine_img, cmap="viridis", transparency=0.5)
add_figure_to_html(
    fig, "Affine registration overlay showing template with affine-transformed subject"
)


# Section 2: Deformable Registration Quality
print("Generating deformable registration comparison...")
fig, axes = plt.subplots(2, 1, figsize=(15, 10))
fig.suptitle(f"Deformable Registration: {subject} to {template}", fontsize=16)

display = plotting.plot_anat(
    template_img, title="Template", display_mode="ortho", figure=fig, axes=axes[0]
)
display = plotting.plot_anat(
    warped_deform_img,
    title="Deformably Warped Subject",
    display_mode="ortho",
    figure=fig,
    axes=axes[1],
)

plt.tight_layout()
add_figure_to_html(
    fig,
    "Deformable registration comparison: template (top) vs. deformably-transformed subject (bottom)",
)

print("Generating deformable overlay...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Deformable Registration Overlay: Template + Subject", fontsize=16)

display = plotting.plot_anat(
    template_img, title="Deformable Overlay", display_mode="z", cut_coords=7, figure=fig
)
display.add_overlay(warped_deform_img, cmap="viridis", transparency=0.5)
add_figure_to_html(
    fig,
    "Deformable registration overlay showing template with deformably-transformed subject",
)


# Section 3: Edge/Contour Overlays
print("Generating affine contours...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Affine Registration with Contours", fontsize=16)

display = plotting.plot_anat(
    warped_affine_img,
    title="Affine with Template Contours",
    display_mode="z",
    cut_coords=7,
    figure=fig,
    dim=-1,
)
display.add_contours(
    template_img, levels=[100, 1000, 2000], colors="r", transparency=0.7
)
add_figure_to_html(
    fig, "Affine-transformed subject with template contours overlaid in red"
)

print("Generating deformable contours...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Deformable Registration with Contours", fontsize=16)

display = plotting.plot_anat(
    warped_deform_img,
    title="Deformable with Template Contours",
    display_mode="z",
    cut_coords=7,
    figure=fig,
    dim=-1,
)
display.add_contours(
    template_img, levels=[100, 1000, 2000], colors="r", transparency=0.7
)
add_figure_to_html(
    fig, "Deformably-transformed subject with template contours overlaid in red"
)


# Section 4: Label Map Visualization
print("Generating atlas segmentation on template...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Atlas Segmentation on Template", fontsize=16)

display = plotting.plot_roi(
    dseg_img,
    bg_img=template_img,
    title="Atlas Labels on Template",
    display_mode="z",
    cut_coords=7,
    cmap="tab20",
    figure=fig,
)
add_figure_to_html(fig, "Atlas segmentation labels overlaid on template image")

print("Generating atlas segmentation on warped subject...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Atlas Segmentation on Warped Subject", fontsize=16)

display = plotting.plot_roi(
    dseg_img,
    bg_img=warped_deform_img,
    title="Atlas Labels on Deformably Warped Subject",
    display_mode="z",
    cut_coords=7,
    cmap="tab20",
    figure=fig,
)
add_figure_to_html(
    fig, "Atlas segmentation labels overlaid on deformably-transformed subject"
)


# Section 5: Warp Field Visualization
print("Calculating warp field magnitude...")
warp_data = warp_img.get_fdata()
warp_stats = {}

if len(warp_data.shape) == 5 and warp_data.shape[-1] == 3:
    # Calculate magnitude of displacement
    warp_magnitude = np.sqrt(np.sum(warp_data**2, axis=-1))

    # Store statistics
    warp_stats = {
        "mean": float(np.mean(warp_magnitude)),
        "std": float(np.std(warp_magnitude)),
        "max": float(np.max(warp_magnitude)),
    }

    print(f"Warp magnitude statistics:")
    print(f"  Mean: {warp_stats['mean']:.3f}")
    print(f"  Std:  {warp_stats['std']:.3f}")
    print(f"  Max:  {warp_stats['max']:.3f}")

    # Create a new NIfTI image for the magnitude
    warp_mag_img = nib.Nifti1Image(warp_magnitude, warp_img.affine, warp_img.header)

    # Visualize warp magnitude
    fig = plt.figure(figsize=(20, 5))
    fig.suptitle(f"Deformation Field Magnitude", fontsize=16)

    display = plotting.plot_stat_map(
        warp_mag_img,
        bg_img=template_img,
        title="Warp Magnitude",
        display_mode="z",
        cut_coords=7,
        cmap="hot",
        figure=fig,
    )
    add_figure_to_html(
        fig,
        f"Deformation field magnitude (mean: {warp_stats['mean']:.2f}, std: {warp_stats['std']:.2f}, max: {warp_stats['max']:.2f})",
    )

    # Visualize individual warp components
    print("Generating warp field components...")
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle(f"Deformation Field Components (X, Y, Z)", fontsize=16)

    component_names = ["X", "Y", "Z"]
    for i, (ax, name) in enumerate(zip(axes, component_names)):
        component_img = nib.Nifti1Image(
            warp_data[..., i], warp_img.affine, warp_img.header
        )
        display = plotting.plot_stat_map(
            component_img,
            bg_img=template_img,
            title=f"Warp {name} component",
            display_mode="z",
            cut_coords=5,
            cmap="coolwarm",
            symmetric_cbar=True,
            axes=ax,
        )
    plt.tight_layout()
    add_figure_to_html(fig, "Deformation field components in X, Y, and Z directions")
else:
    print(f"Warp field has unexpected shape: {warp_data.shape}")


# Generate HTML Report
print("Generating HTML report...")

html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Registration QC Report - {subject}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #3498db;
            padding-left: 10px;
        }}
        .info-box {{
            background-color: #ecf0f1;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 20px 0;
        }}
        .info-box p {{
            margin: 5px 0;
        }}
        .figure {{
            margin: 30px 0;
            text-align: center;
        }}
        .figure img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #ddd;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .caption {{
            margin-top: 10px;
            font-style: italic;
            color: #666;
        }}
        .stats-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        .stats-table th, .stats-table td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        .stats-table th {{
            background-color: #3498db;
            color: white;
        }}
        .stats-table tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        .checklist {{
            list-style-type: none;
            padding-left: 0;
        }}
        .checklist li::before {{
            content: "✓ ";
            color: #27ae60;
            font-weight: bold;
            margin-right: 5px;
        }}
        .summary-box {{
            background-color: #e8f8f5;
            border: 2px solid #27ae60;
            border-radius: 5px;
            padding: 20px;
            margin: 30px 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Registration Quality Control Report</h1>
        
        <div class="info-box">
            <p><strong>Subject ID:</strong> {subject}</p>
            <p><strong>Template:</strong> {template}</p>
            <p><strong>Stain:</strong> {stain}</p>
            <p><strong>Generated:</strong> {Path(output_html).name}</p>
        </div>

        <h2>Image Information</h2>
        <table class="stats-table">
            <tr>
                <th>Image Type</th>
                <th>Shape</th>
            </tr>
            <tr>
                <td>Template</td>
                <td>{image_info['template_shape']}</td>
            </tr>
            <tr>
                <td>Subject (original)</td>
                <td>{image_info['subject_shape']}</td>
            </tr>
            <tr>
                <td>Affine Warped</td>
                <td>{image_info['warped_affine_shape']}</td>
            </tr>
            <tr>
                <td>Deformably Warped</td>
                <td>{image_info['warped_deform_shape']}</td>
            </tr>
            <tr>
                <td>Warp Field</td>
                <td>{image_info['warp_shape']}</td>
            </tr>
            <tr>
                <td>Atlas Segmentation</td>
                <td>{image_info['dseg_shape']}</td>
            </tr>
        </table>
"""

# Add warp statistics if available
if warp_stats:
    html_content += f"""
        <h2>Deformation Field Statistics</h2>
        <table class="stats-table">
            <tr>
                <th>Statistic</th>
                <th>Value</th>
            </tr>
            <tr>
                <td>Mean Displacement</td>
                <td>{warp_stats['mean']:.3f}</td>
            </tr>
            <tr>
                <td>Standard Deviation</td>
                <td>{warp_stats['std']:.3f}</td>
            </tr>
            <tr>
                <td>Maximum Displacement</td>
                <td>{warp_stats['max']:.3f}</td>
            </tr>
        </table>
"""

# Add section descriptions
sections = [
    (
        "1. Affine Registration Quality",
        "Comparison of the template and affine-transformed subject images. Affine registration accounts for global transformations including translation, rotation, scaling, and shearing.",
    ),
    (
        "2. Deformable Registration Quality",
        "Comparison of the template and deformably-transformed subject images. Deformable registration accounts for local anatomical variations through non-linear warping.",
    ),
    (
        "3. Edge/Contour Overlays",
        "Visualization of registration quality using edge contours. The red contours from the template are overlaid on the transformed subject images to assess alignment accuracy.",
    ),
    (
        "4. Label Map Visualization",
        "Atlas segmentation labels overlaid on both template and transformed subject images. This shows how well anatomical regions align after registration.",
    ),
    (
        "5. Warp Field Visualization",
        "Visualization of the deformation field showing the magnitude and directional components of the transformation applied during deformable registration.",
    ),
]

section_idx = 0
for i, fig_data in enumerate(html_figures):
    # Add section header before the first figure of each section
    if i in [0, 2, 4, 6, 8, 10]:  # Indices where new sections start
        if section_idx < len(sections):
            html_content += f"""
        <h2>{sections[section_idx][0]}</h2>
        <p>{sections[section_idx][1]}</p>
"""
            section_idx += 1

    html_content += f"""
        <div class="figure">
            <img src="{fig_data['image']}" alt="Registration QC Figure {i+1}">
            <div class="caption">{fig_data['caption']}</div>
        </div>
"""

# Add summary
html_content += f"""
        <div class="summary-box">
            <h2>Summary</h2>
            <p>Registration QC report completed successfully for subject <strong>{subject}</strong> registered to template <strong>{template}</strong> using stain <strong>{stain}</strong>.</p>
            
            <h3>QC Visualizations Generated:</h3>
            <ul class="checklist">
                <li>Affine registration comparison</li>
                <li>Affine registration overlay</li>
                <li>Deformable registration comparison</li>
                <li>Deformable registration overlay</li>
                <li>Affine registration with contours</li>
                <li>Deformable registration with contours</li>
                <li>Atlas segmentation on template</li>
                <li>Atlas segmentation on warped subject</li>
                <li>Deformation field magnitude</li>
                <li>Deformation field components</li>
            </ul>
        </div>
    </div>
</body>
</html>
"""

# Write HTML to output file
with open(output_html, "w") as f:
    f.write(html_content)

print("=" * 80)
print(f"Registration QC Report Summary")
print("=" * 80)
print(f"Subject: {subject}")
print(f"Template: {template}")
print(f"Stain: {stain}")
print(f"Output HTML: {output_html}")
print("=" * 80)
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
print("=" * 80)
print(f"HTML report written to: {output_html}")
