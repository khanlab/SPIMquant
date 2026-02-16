#!/usr/bin/env python3
"""
Generate modular Mermaid DAG diagrams from SPIMquant workflow.

This script:
1. Generates the full rulegraph from Snakemake
2. Splits it into functional subgraphs based on workflow stages
3. Creates individual mermaid files for each stage
4. Renders diagrams to SVG format

Usage:
    python generate_dag_diagrams.py [--bids-dir PATH] [--output-dir PATH]
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple


# Define workflow stages based on rule naming patterns and dependencies
WORKFLOW_STAGES = {
    "01_import": {
        "description": "Import and setup templates, masks, and atlases",
        "patterns": [
            r"^import_",
            r"^generic_lut_",
            r"^copy_template_",
        ],
    },
    "02_preprocessing": {
        "description": "Image preprocessing and downsampling",
        "patterns": [
            r"^get_downsampled_",
            r"^pre_atropos$",
        ],
    },
    "03_masking": {
        "description": "Brain masking using atropos segmentation",
        "patterns": [
            r"^atropos_",
            r"^post_atropos$",
            r"^create_mask_",
            r"^affine_transform_template_mask",
        ],
    },
    "04_correction": {
        "description": "Intensity correction and normalization",
        "patterns": [
            r"^n4$",
            r"^n4_biasfield$",
            r"^apply_mask_to_corrected$",
        ],
    },
    "05_registration": {
        "description": "Template registration (affine and deformable)",
        "patterns": [
            r"^init_affine_",
            r"^affine_reg$",
            r"^deform_reg$",
            r"^convert_ras_to_itk$",
        ],
    },
    "06_transform": {
        "description": "Apply transformations to images and atlases",
        "patterns": [
            r"^deform_spim_nii_to_template",
            r"^deform_template_dseg_to_subject",
            r"^transform_regionprops",
        ],
    },
    "07_segmentation": {
        "description": "Segmentation of pathology (threshold, multi-otsu)",
        "patterns": [
            r"^threshold$",
            r"^multiotsu$",
            r"^fieldfrac$",
            r"^deform_fieldfrac",
        ],
    },
    "08_quantification": {
        "description": "Region properties and quantification",
        "patterns": [
            r"^compute_filtered_regionprops$",
            r"^map_regionprops_to_atlas",
            r"^counts_per_voxel",
            r"^aggregate_regionprops",
        ],
    },
    "09_statistics": {
        "description": "Statistical analysis and atlas mapping",
        "patterns": [
            r"^merge_.*segstats",
            r"^map_.*_to_.*_nii$",
            r"^map_img_to_roi",
        ],
    },
    "10_qc": {
        "description": "Quality control and reporting",
        "patterns": [
            r"^registration_qc_report$",
        ],
    },
    "11_patches": {
        "description": "Extract image patches for analysis",
        "patterns": [
            r"^create_.*_patches$",
        ],
    },
}


def run_command(cmd: List[str], capture_output=True) -> Tuple[int, str, str]:
    """Run a command and return exit code, stdout, stderr."""
    result = subprocess.run(
        cmd, capture_output=capture_output, text=True, cwd=Path(__file__).parent.parent.parent
    )
    return result.returncode, result.stdout, result.stderr


def generate_rulegraph(bids_dir: Path, output_dir: Path, template: str = "ABAv3") -> str:
    """Generate the full rulegraph using Snakemake."""
    print(f"Generating rulegraph from {bids_dir}...")
    
    cmd = [
        sys.executable,
        "-m",
        "spimquant.run",
        str(bids_dir),
        str(output_dir),
        "participant",
        "--sloppy",
        "--template",
        template,
        "--rulegraph",
        "mermaid-js",
    ]
    
    returncode, stdout, stderr = run_command(cmd)
    
    if returncode != 0:
        print(f"Error generating rulegraph: {stderr}", file=sys.stderr)
        sys.exit(1)
    
    # Extract just the mermaid graph (skip app messages)
    lines = stdout.split("\n")
    mermaid_start = None
    for i, line in enumerate(lines):
        if line.strip() == "---" and i + 1 < len(lines) and lines[i + 1].strip().startswith("title:"):
            mermaid_start = i
            break
    
    if mermaid_start is None:
        print("Could not find mermaid graph in output", file=sys.stderr)
        sys.exit(1)
    
    return "\n".join(lines[mermaid_start:])


def parse_mermaid_graph(mermaid: str) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """
    Parse mermaid flowchart to extract nodes and edges.
    
    Returns:
        nodes: Dict mapping node id to node label
        edges: List of (source_id, target_id) tuples
    """
    nodes = {}
    edges = []
    
    for line in mermaid.split("\n"):
        line = line.strip()
        
        # Parse node definition: id0[label]
        node_match = re.match(r"(id\d+)\[(.*?)\]", line)
        if node_match:
            node_id = node_match.group(1)
            node_label = node_match.group(2)
            nodes[node_id] = node_label
            continue
        
        # Parse edge: id1 --> id2
        edge_match = re.match(r"(id\d+)\s+-->\s+(id\d+)", line)
        if edge_match:
            source = edge_match.group(1)
            target = edge_match.group(2)
            edges.append((source, target))
    
    return nodes, edges


def classify_node(node_label: str, stages: Dict) -> str:
    """Classify a node into a workflow stage based on patterns."""
    for stage_name, stage_info in stages.items():
        for pattern in stage_info["patterns"]:
            if re.search(pattern, node_label):
                return stage_name
    return "other"


def get_node_dependencies(node_id: str, edges: List[Tuple[str, str]]) -> Set[str]:
    """Get all nodes that this node depends on (directly or indirectly)."""
    dependencies = set()
    to_process = {node_id}
    
    while to_process:
        current = to_process.pop()
        for source, target in edges:
            if target == current and source not in dependencies:
                dependencies.add(source)
                to_process.add(source)
    
    return dependencies


def create_subgraph(
    stage_name: str,
    stage_info: Dict,
    nodes: Dict[str, str],
    edges: List[Tuple[str, str]],
    all_classifications: Dict[str, str],
) -> str:
    """Create a mermaid subgraph for a specific workflow stage."""
    # Find nodes in this stage
    stage_nodes = {
        node_id: label
        for node_id, label in nodes.items()
        if all_classifications[node_id] == stage_name
    }
    
    if not stage_nodes:
        return None
    
    # Find edges within this stage or connecting to it
    stage_node_ids = set(stage_nodes.keys())
    relevant_edges = []
    
    for source, target in edges:
        # Include edges where at least one node is in this stage
        if source in stage_node_ids or target in stage_node_ids:
            relevant_edges.append((source, target))
    
    # Build mermaid graph
    lines = [
        "---",
        f"title: {stage_info['description']}",
        "---",
        "flowchart TB",
    ]
    
    # Add all relevant nodes (including dependencies from other stages)
    all_relevant_nodes = set()
    for source, target in relevant_edges:
        all_relevant_nodes.add(source)
        all_relevant_nodes.add(target)
    
    for node_id in sorted(all_relevant_nodes):
        if node_id in nodes:
            label = nodes[node_id]
            lines.append(f"\t{node_id}[{label}]")
    
    # Add styling for stage nodes vs dependency nodes
    for node_id in sorted(all_relevant_nodes):
        if node_id in stage_node_ids:
            # Primary stage nodes - bold style
            lines.append(f"\tstyle {node_id} fill:#57D996,stroke:#2D8659,stroke-width:3px")
        else:
            # Dependency nodes from other stages - lighter style
            lines.append(f"\tstyle {node_id} fill:#E8F4F8,stroke:#AAA,stroke-width:1px")
    
    # Add edges
    for source, target in relevant_edges:
        lines.append(f"\t{source} --> {target}")
    
    return "\n".join(lines)


def render_mermaid_to_svg(mermaid_file: Path, output_file: Path):
    """Render a mermaid file to SVG using mmdc CLI if available."""
    try:
        # Check if mmdc (mermaid-cli) is available
        result = subprocess.run(
            ["mmdc", "--version"], capture_output=True, text=True
        )
        if result.returncode == 0:
            print(f"  Rendering {mermaid_file.name} to SVG...")
            subprocess.run(
                ["mmdc", "-i", str(mermaid_file), "-o", str(output_file), "-t", "neutral"],
                check=True,
            )
            return True
    except FileNotFoundError:
        pass
    
    return False


def main():
    parser = argparse.ArgumentParser(
        description="Generate modular Mermaid DAG diagrams for SPIMquant workflow"
    )
    parser.add_argument(
        "--bids-dir",
        type=Path,
        help="Path to BIDS dataset (default: tests/bids_ds)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Path to output directory for workflow (default: /tmp/spimquant_output)",
    )
    parser.add_argument(
        "--template",
        type=str,
        default="ABAv3",
        help="Template to use for registration (default: ABAv3)",
    )
    args = parser.parse_args()
    
    # Get repository root
    repo_root = Path(__file__).parent.parent.parent
    
    # Set defaults
    bids_dir = args.bids_dir or repo_root / "tests" / "bids_ds"
    output_dir = args.output_dir or Path("/tmp/spimquant_output")
    
    # Create output directories
    figures_dir = repo_root / "docs" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("SPIMquant DAG Diagram Generator")
    print("=" * 80)
    print(f"BIDS directory: {bids_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Figures directory: {figures_dir}")
    print()
    
    # Generate full rulegraph
    full_rulegraph = generate_rulegraph(bids_dir, output_dir, args.template)
    
    # Save full rulegraph
    full_rulegraph_file = figures_dir / "rulegraph_full.mermaid"
    print(f"Saving full rulegraph to {full_rulegraph_file}...")
    full_rulegraph_file.write_text(full_rulegraph)
    
    # Parse the graph
    print("Parsing rulegraph...")
    nodes, edges = parse_mermaid_graph(full_rulegraph)
    print(f"  Found {len(nodes)} nodes and {len(edges)} edges")
    
    # Classify all nodes
    print("Classifying nodes into workflow stages...")
    classifications = {}
    for node_id, label in nodes.items():
        stage = classify_node(label, WORKFLOW_STAGES)
        classifications[node_id] = stage
        if stage == "other":
            print(f"  Warning: Node '{label}' not classified into any stage")
    
    # Print classification summary
    print("\nClassification summary:")
    for stage_name, stage_info in WORKFLOW_STAGES.items():
        count = sum(1 for c in classifications.values() if c == stage_name)
        if count > 0:
            print(f"  {stage_name}: {count} nodes - {stage_info['description']}")
    
    # Generate subgraphs for each stage
    print("\nGenerating stage-specific diagrams...")
    rendered_count = 0
    for stage_name, stage_info in WORKFLOW_STAGES.items():
        print(f"\n{stage_name}: {stage_info['description']}")
        
        subgraph = create_subgraph(
            stage_name, stage_info, nodes, edges, classifications
        )
        
        if subgraph:
            # Save mermaid file
            mermaid_file = figures_dir / f"dag_{stage_name}.mermaid"
            mermaid_file.write_text(subgraph)
            print(f"  Saved: {mermaid_file.name}")
            
            # Try to render to SVG
            svg_file = figures_dir / f"dag_{stage_name}.svg"
            if render_mermaid_to_svg(mermaid_file, svg_file):
                print(f"  Rendered: {svg_file.name}")
                rendered_count += 1
        else:
            print(f"  Skipped (no nodes in this stage)")
    
    print("\n" + "=" * 80)
    print("Summary:")
    print(f"  Full rulegraph: {full_rulegraph_file}")
    print(f"  Stage diagrams: {len([s for s in WORKFLOW_STAGES if (figures_dir / f'dag_{s}.mermaid').exists()])} created")
    
    if rendered_count == 0:
        print("\nNote: To render diagrams to SVG/PNG, install mermaid-cli:")
        print("  npm install -g @mermaid-js/mermaid-cli")
    else:
        print(f"  SVG renders: {rendered_count} created")
    
    print("=" * 80)


if __name__ == "__main__":
    main()
