"""
MkDocs hooks for SPIMquant documentation.

This module provides hooks that run during the MkDocs build process.
"""

import subprocess
import sys
from pathlib import Path


def on_pre_build(config):
    """
    Hook that runs before the build starts.

    Generates CLI reference documentation from argparse.
    """
    print("Generating CLI reference documentation...")

    # Path to the gen_cli_docs script
    docs_root = Path(__file__).parent
    script_path = docs_root / "scripts" / "gen_cli_docs.py"

    try:
        # Run the script to generate CLI docs
        result = subprocess.run(
            [sys.executable, str(script_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        print("CLI reference documentation generated successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error generating CLI documentation: {e}", file=sys.stderr)
        print(f"stdout: {e.stdout}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        # Don't fail the build, just warn
        print("Warning: CLI documentation generation failed. Continuing build...")
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        print("Warning: CLI documentation generation failed. Continuing build...")
