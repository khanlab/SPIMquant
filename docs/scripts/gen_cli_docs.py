#!/usr/bin/env python3
"""
Generate CLI documentation from argparse parser.

This script generates markdown documentation for the SPIMquant CLI
by extracting information from the argparse parser defined in spimquant/run.py.
"""

import argparse
import sys
from pathlib import Path


def format_choices(choices):
    """Format choices for display."""
    if choices is None:
        return ""
    return f" (choices: {', '.join(map(str, choices))})"


def format_default(default):
    """Format default value for display."""
    if default is None or default == argparse.SUPPRESS:
        return ""
    if isinstance(default, list):
        if not default:
            return ""
        return f" (default: {', '.join(map(str, default))})"
    return f" (default: {default})"


def format_action_type(action):
    """Get a human-readable type for the action."""
    if action.type is not None:
        return f" ({action.type.__name__})"
    elif isinstance(action, argparse._StoreAction):
        return " (str)"
    return ""


def generate_cli_docs(parser, output_file):
    """
    Generate markdown documentation from argparse parser.

    Args:
        parser: ArgumentParser instance
        output_file: Path to output markdown file
    """
    output = []
    output.append("# Command Line Interface Reference\n")
    output.append(
        "This page provides auto-generated documentation for all SPIMquant command-line options.\n"
    )
    output.append("## Overview\n")
    output.append("```bash")
    output.append("spimquant <bids_dir> <output_dir> <analysis_level> [options]")
    output.append("```\n")

    # Get description
    if parser.description:
        output.append("## Description\n")
        output.append(f"{parser.description}\n")

    # Positional arguments
    positionals = []
    optionals = []

    for action in parser._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        if isinstance(action, argparse._VersionAction):
            continue

        if action.option_strings:
            optionals.append(action)
        else:
            positionals.append(action)

    # Document positional arguments
    if positionals:
        output.append("## Positional Arguments\n")
        for action in positionals:
            if action.dest == "==SUPPRESS==":
                continue
            output.append(f"### `{action.dest}`\n")
            if action.help:
                help_text = action.help.replace("\n", " ").strip()
                output.append(f"{help_text}")
            output.append(f"{format_action_type(action)}")
            output.append(f"{format_choices(action.choices)}")
            output.append("\n")

    # Document optional arguments
    if optionals:
        output.append("## Optional Arguments\n")

        # Group by category if possible (based on argument group)
        # For now, list all optional arguments
        for action in optionals:
            if action.dest == "==SUPPRESS==":
                continue

            # Format option strings
            opts = ", ".join(f"`{opt}`" for opt in action.option_strings)
            output.append(f"### {opts}\n")

            if action.help:
                help_text = action.help.replace("%(default)s", str(action.default))
                help_text = help_text.replace("\n", " ").strip()
                output.append(f"{help_text}")

            type_info = format_action_type(action)
            if type_info:
                output.append(type_info)

            choices = format_choices(action.choices)
            if choices:
                output.append(choices)

            default = format_default(action.default)
            if default:
                output.append(default)

            # Handle special actions
            if isinstance(action, argparse._StoreTrueAction):
                output.append(" (flag)")
            elif isinstance(action, argparse._StoreFalseAction):
                output.append(" (flag)")

            if action.nargs:
                if action.nargs == "+":
                    output.append(" (accepts one or more values)")
                elif action.nargs == "*":
                    output.append(" (accepts zero or more values)")
                elif isinstance(action.nargs, int):
                    output.append(f" (accepts {action.nargs} values)")

            output.append("\n")

    # Write to file
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(output))
    print(f"Generated CLI documentation at {output_path}")


def main():
    """Main entry point."""
    # Import parser from spimquant
    try:
        # Add parent directory to path
        project_root = Path(__file__).parent.parent.parent
        sys.path.insert(0, str(project_root))

        from spimquant.run import get_parser

        parser = get_parser()

        # Default output location
        output_file = Path(__file__).parent.parent / "reference" / "cli_reference.md"

        # Allow override from command line
        if len(sys.argv) > 1:
            output_file = Path(sys.argv[1])

        generate_cli_docs(parser, output_file)

    except ImportError as e:
        print(f"Error importing parser: {e}", file=sys.stderr)
        print("Make sure all dependencies are installed.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error generating documentation: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
