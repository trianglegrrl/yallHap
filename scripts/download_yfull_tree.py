#!/usr/bin/env python3
"""
Download YFull tree from GitHub.

This script downloads the latest YFull tree JSON file
from the official YFullTeam GitHub repository.
"""

import argparse
import sys
from pathlib import Path

import requests


YFULL_TREE_URL = "https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download YFull tree from GitHub"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("data/yfull_tree.json"),
        help="Output file path [default: data/yfull_tree.json]"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="Request timeout in seconds [default: 60]"
    )

    args = parser.parse_args()

    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading YFull tree from {YFULL_TREE_URL}...")

    try:
        response = requests.get(YFULL_TREE_URL, timeout=args.timeout)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"Error downloading tree: {e}", file=sys.stderr)
        return 1

    # Validate JSON
    try:
        import json
        data = response.json()
        print(f"  Tree contains {len(data)} nodes")
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON: {e}", file=sys.stderr)
        return 1

    # Write to file
    with open(args.output, "wb") as f:
        f.write(response.content)

    print(f"Saved to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
