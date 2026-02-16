"""
MRI to SPIM Registration Quality Control Report Generator

This script generates an HTML report with comprehensive quality control 
visualizations for MRI to SPIM registration, including comparisons of SPIM 
and transformed MRI images, overlay visualizations, and warp field visualizations.
"""

import base64
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from nilearn import plotting

# Load inputs from snakemake object
spim_path = snakemake.input.spim
mri_path = snakemake.input.mri
warped_affine_path = snakemake.input.warped_affine
warped_deform_path = snakemake.input.warped_deform
warp_field_path = snakemake.input.warp

output_html = snakemake.output.report_html

# Get wildcards for report title
subject = snakemake.wildcards.subject
stain = snakemake.wildcards.stain

print(
    f"Processing MRI to SPIM registration QC for subject {subject}, stain {stain}"
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
spim_img = nib.load(spim_path)
mri_img = nib.load(mri_path)
warped_affine_img = nib.load(warped_affine_path)
warped_deform_img = nib.load(warped_deform_path)
warp_img = nib.load(warp_field_path)

print(f"SPIM shape: {spim_img.shape}")
print(f"MRI shape: {mri_img.shape}")
print(f"Warped (affine) shape: {warped_affine_img.shape}")
print(f"Warped (deform) shape: {warped_deform_img.shape}")
print(f"Warp field shape: {warp_img.shape}")

# Store image info for report
image_info = {
    "spim_shape": spim_img.shape,
    "mri_shape": mri_img.shape,
    "warped_affine_shape": warped_affine_img.shape,
    "warped_deform_shape": warped_deform_img.shape,
    "warp_shape": warp_img.shape,
}

# Section 1: Affine Registration Quality
print("Generating affine registration comparison...")
fig, axes = plt.subplots(2, 1, figsize=(15, 10))
fig.suptitle(f"Affine Registration: MRI to SPIM (subject {subject})", fontsize=16)

display = plotting.plot_anat(
    spim_img,
    title="SPIM Reference",
    display_mode="ortho",
    figure=fig,
    axes=axes[0],
    cut_coords=(0, 0, 0),
)
display = plotting.plot_anat(
    warped_affine_img,
    title="Affine Warped MRI",
    display_mode="ortho",
    figure=fig,
    axes=axes[1],
    cut_coords=(0, 0, 0),
    dim=-1.5,
)

plt.tight_layout()
add_figure_to_html(
    fig,
    "Affine registration comparison: SPIM reference (top) vs. affine-transformed MRI (bottom)",
)

print("Generating affine overlay...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Affine Registration Overlay: SPIM + MRI", fontsize=16)

display = plotting.plot_anat(
    spim_img, title="Affine Overlay", display_mode="z", cut_coords=7, figure=fig
)
display.add_overlay(warped_affine_img, cmap="viridis", transparency=0.5)
add_figure_to_html(
    fig, "Affine registration overlay showing SPIM with affine-transformed MRI"
)


# Section 2: Deformable Registration Quality
print("Generating deformable registration comparison...")
fig, axes = plt.subplots(2, 1, figsize=(15, 10))
fig.suptitle(f"Deformable Registration: MRI to SPIM (subject {subject})", fontsize=16)

display = plotting.plot_anat(
    spim_img,
    title="SPIM Reference",
    display_mode="ortho",
    cut_coords=(0, 0, 0),
    figure=fig,
    axes=axes[0],
)
display = plotting.plot_anat(
    warped_deform_img,
    title="Deformably Warped MRI",
    display_mode="ortho",
    cut_coords=(0, 0, 0),
    dim=-1.5,
    figure=fig,
    axes=axes[1],
)

plt.tight_layout()
add_figure_to_html(
    fig,
    "Deformable registration comparison: SPIM reference (top) vs. deformably-transformed MRI (bottom)",
)

print("Generating deformable overlay...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Deformable Registration Overlay: SPIM + MRI", fontsize=16)

display = plotting.plot_anat(
    spim_img, title="Deformable Overlay", display_mode="z", cut_coords=7, figure=fig
)
display.add_overlay(warped_deform_img, cmap="viridis", transparency=0.5)
add_figure_to_html(
    fig,
    "Deformable registration overlay showing SPIM with deformably-transformed MRI",
)


# Section 3: Edge/Contour Overlays
print("Generating affine contours...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Affine Registration with Contours", fontsize=16)

display = plotting.plot_anat(
    warped_affine_img,
    title="Affine with SPIM Contours",
    display_mode="z",
    cut_coords=7,
    figure=fig,
    dim=-1,
)
display.add_contours(
    spim_img, levels=[100, 1000, 2000], colors="r", transparency=0.7
)
add_figure_to_html(
    fig, "Affine-transformed MRI with SPIM contours overlaid in red"
)

print("Generating deformable contours...")
fig = plt.figure(figsize=(20, 5))
fig.suptitle(f"Deformable Registration with Contours", fontsize=16)

display = plotting.plot_anat(
    warped_deform_img,
    title="Deformable with SPIM Contours",
    display_mode="z",
    cut_coords=7,
    figure=fig,
    dim=-1,
)
display.add_contours(
    spim_img, levels=[100, 1000, 2000], colors="r", transparency=0.7
)
add_figure_to_html(
    fig, "Deformably-transformed MRI with SPIM contours overlaid in red"
)


# Section 4: Warp Field Visualization
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
    fig.suptitle(f"Deformation Field Magnitude (MRI to SPIM)", fontsize=16)

    display = plotting.plot_stat_map(
        warp_mag_img,
        bg_img=spim_img,
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

    # Visualize individual displacement components
    print("Generating displacement component visualizations...")
    fig, axes = plt.subplots(1, 3, figsize=(20, 5))
    fig.suptitle(f"Deformation Field Components (MRI to SPIM)", fontsize=16)

    component_names = ["X (Left-Right)", "Y (Anterior-Posterior)", "Z (Superior-Inferior)"]
    
    for i, (ax, comp_name) in enumerate(zip(axes, component_names)):
        # Create NIfTI image for each component
        comp_img = nib.Nifti1Image(
            warp_data[..., 0, i], warp_img.affine, warp_img.header
        )
        
        display = plotting.plot_stat_map(
            comp_img,
            bg_img=spim_img,
            title=comp_name,
            display_mode="z",
            cut_coords=[0],
            cmap="RdBu_r",
            axes=ax,
        )
    
    plt.tight_layout()
    add_figure_to_html(
        fig, "Directional components of the deformation field showing displacement in each axis"
    )

else:
    print(f"Warp field has unexpected shape: {warp_data.shape}")


# Generate HTML Report
print("Generating HTML report...")

html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MRI to SPIM Registration QC Report - {subject}</title>
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
        <h1>MRI to SPIM Registration Quality Control Report</h1>
        
        <div class="info-box">
            <p><strong>Subject ID:</strong> {subject}</p>
            <p><strong>Reference Stain:</strong> {stain}</p>
            <p><strong>Registration Type:</strong> MRI to SPIM</p>
            <p><strong>Generated:</strong> {Path(output_html).name}</p>
        </div>

        <h2>Image Information</h2>
        <table class="stats-table">
            <tr>
                <th>Image Type</th>
                <th>Shape</th>
            </tr>
            <tr>
                <td>SPIM Reference</td>
                <td>{image_info['spim_shape']}</td>
            </tr>
            <tr>
                <td>MRI (original)</td>
                <td>{image_info['mri_shape']}</td>
            </tr>
            <tr>
                <td>Affine Warped MRI</td>
                <td>{image_info['warped_affine_shape']}</td>
            </tr>
            <tr>
                <td>Deformably Warped MRI</td>
                <td>{image_info['warped_deform_shape']}</td>
            </tr>
            <tr>
                <td>Warp Field</td>
                <td>{image_info['warp_shape']}</td>
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
        "Comparison of the SPIM reference and affine-transformed MRI images. Affine registration accounts for global transformations including translation, rotation, scaling, and shearing.",
    ),
    (
        "2. Deformable Registration Quality",
        "Comparison of the SPIM reference and deformably-transformed MRI images. Deformable registration accounts for local anatomical variations through non-linear warping.",
    ),
    (
        "3. Edge/Contour Overlays",
        "Visualization of registration quality using edge contours. The red contours from the SPIM reference are overlaid on the transformed MRI images to assess alignment accuracy.",
    ),
    (
        "4. Warp Field Visualization",
        "Visualization of the deformation field showing the magnitude and directional components of the transformation applied during deformable registration from MRI to SPIM space.",
    ),
]

section_idx = 0
for i, fig_data in enumerate(html_figures):
    # Add section header before the first figure of each section
    if i in [0, 2, 4, 6]:  # Indices where new sections start
        if section_idx < len(sections):
            html_content += f"""
        <h2>{sections[section_idx][0]}</h2>
        <p>{sections[section_idx][1]}</p>
"""
            section_idx += 1

    html_content += f"""
        <div class="figure">
            <img src="{fig_data['image']}" alt="MRI to SPIM Registration QC Figure {i+1}">
            <div class="caption">{fig_data['caption']}</div>
        </div>
"""

# Add summary
html_content += f"""
        <div class="summary-box">
            <h2>Summary</h2>
            <p>MRI to SPIM registration QC report completed successfully for subject <strong>{subject}</strong> using SPIM stain <strong>{stain}</strong> as reference.</p>
            
            <h3>QC Visualizations Generated:</h3>
            <ul class="checklist">
                <li>Affine registration comparison</li>
                <li>Affine registration overlay</li>
                <li>Deformable registration comparison</li>
                <li>Deformable registration overlay</li>
                <li>Affine registration with contours</li>
                <li>Deformable registration with contours</li>
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
print(f"MRI to SPIM Registration QC Report Summary")
print("=" * 80)
print(f"Subject: {subject}")
print(f"Reference Stain: {stain}")
print(f"Output HTML: {output_html}")
print("=" * 80)
print("QC visualizations generated:")
print("  ✓ Affine registration comparison")
print("  ✓ Affine registration overlay")
print("  ✓ Deformable registration comparison")
print("  ✓ Deformable registration overlay")
print("  ✓ Affine registration with contours")
print("  ✓ Deformable registration with contours")
print("  ✓ Deformation field magnitude")
print("  ✓ Deformation field components")
print("=" * 80)
print(f"HTML report written to: {output_html}")
