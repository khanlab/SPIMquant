#!/usr/bin/env python3
"""
Render Mermaid diagrams to SVG/PNG format.

This script renders all .mermaid files in docs/figures/ to SVG and optionally PNG format.
It handles the puppeteer sandbox configuration needed for CI environments.

Usage:
    python3 render_diagrams.py [--format svg|png|both] [--theme neutral]
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def check_mmdc_available():
    """Check if mermaid-cli (mmdc) is installed."""
    try:
        result = subprocess.run(
            ["mmdc", "--version"], capture_output=True, text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def render_diagram(mermaid_file: Path, output_format: str = "svg", theme: str = "neutral"):
    """
    Render a single mermaid diagram to the specified format.
    
    Args:
        mermaid_file: Path to .mermaid file
        output_format: Output format ('svg' or 'png')
        theme: Mermaid theme to use ('neutral', 'default', 'dark', 'forest', etc.)
    
    Returns:
        True if successful, False otherwise
    """
    output_file = mermaid_file.with_suffix(f".{output_format}")
    
    # Create puppeteer config for CI environments (no sandbox)
    puppeteer_config = json.dumps({
        "args": ["--no-sandbox", "--disable-setuid-sandbox"]
    })
    
    try:
        result = subprocess.run(
            [
                "mmdc",
                "-i", str(mermaid_file),
                "-o", str(output_file),
                "-t", theme,
                "--puppeteerConfigFile", "/dev/stdin",
            ],
            input=puppeteer_config,
            capture_output=True,
            text=True,
        )
        
        if result.returncode == 0:
            return True
        else:
            print(f"  Error rendering {mermaid_file.name}: {result.stderr}", file=sys.stderr)
            return False
            
    except Exception as e:
        print(f"  Exception rendering {mermaid_file.name}: {e}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Render Mermaid diagrams to SVG/PNG format"
    )
    parser.add_argument(
        "--format",
        choices=["svg", "png", "both"],
        default="svg",
        help="Output format (default: svg)",
    )
    parser.add_argument(
        "--theme",
        default="neutral",
        help="Mermaid theme to use for rendering (default: neutral). The 'neutral' theme provides clean, professional diagrams suitable for documentation.",
    )
    args = parser.parse_args()
    
    # Check if mmdc is available
    if not check_mmdc_available():
        print("Error: mermaid-cli (mmdc) is not installed.", file=sys.stderr)
        print("Install with: npm install -g @mermaid-js/mermaid-cli", file=sys.stderr)
        sys.exit(1)
    
    # Get figures directory
    repo_root = Path(__file__).parent.parent.parent
    figures_dir = repo_root / "docs" / "figures"
    
    if not figures_dir.exists():
        print(f"Error: Figures directory not found: {figures_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Find all mermaid files
    mermaid_files = sorted(figures_dir.glob("*.mermaid"))
    
    if not mermaid_files:
        print(f"No .mermaid files found in {figures_dir}", file=sys.stderr)
        sys.exit(1)
    
    print("=" * 80)
    print("Mermaid Diagram Renderer")
    print("=" * 80)
    print(f"Figures directory: {figures_dir}")
    print(f"Output format: {args.format}")
    print(f"Theme: {args.theme}")
    print(f"Found {len(mermaid_files)} diagram(s) to render")
    print()
    
    # Determine formats to render
    formats = ["svg", "png"] if args.format == "both" else [args.format]
    
    # Render each diagram
    success_count = 0
    fail_count = 0
    
    for mermaid_file in mermaid_files:
        for fmt in formats:
            print(f"Rendering {mermaid_file.name} to {fmt.upper()}...", end=" ")
            if render_diagram(mermaid_file, fmt, args.theme):
                print("✓")
                success_count += 1
            else:
                print("✗")
                fail_count += 1
    
    print()
    print("=" * 80)
    print(f"Rendering complete: {success_count} successful, {fail_count} failed")
    print("=" * 80)
    
    if fail_count > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
