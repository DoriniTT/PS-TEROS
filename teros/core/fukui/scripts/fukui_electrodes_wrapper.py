#!/usr/bin/env python
"""Wrapper script for FukuiGrid.fukui_electrodes() function.

This script is called by the run_fukui_electrodes_calcfunc to compute
the Fukui potential using the electrodes method.

Usage:
    python fukui_electrodes_wrapper.py <chgcar_neutral> <chgcar_fukui> <epsilon>

Arguments:
    chgcar_neutral: Path to CHGCAR file of neutral slab (e.g., CHGCAR_0.00)
    chgcar_fukui: Path to Fukui function file (e.g., CHGCAR_FUKUI.vasp)
    epsilon: Dielectric constant value (float)

Output:
    LOCPOT_FUKUI.vasp in the current working directory
"""

import sys
from pathlib import Path


def main():
    """Run FukuiGrid fukui_electrodes with command-line arguments."""
    if len(sys.argv) != 4:
        print(
            "Usage: fukui_electrodes_wrapper.py <chgcar_neutral> <chgcar_fukui> <epsilon>",
            file=sys.stderr,
        )
        print("\nArguments:", file=sys.stderr)
        print("  chgcar_neutral: Path to CHGCAR of neutral slab", file=sys.stderr)
        print("  chgcar_fukui: Path to Fukui function file", file=sys.stderr)
        print("  epsilon: Dielectric constant value", file=sys.stderr)
        sys.exit(1)

    chgcar_neutral = sys.argv[1]
    chgcar_fukui = sys.argv[2]
    try:
        epsilon = float(sys.argv[3])
    except ValueError:
        print(f"Error: epsilon must be a number, got '{sys.argv[3]}'", file=sys.stderr)
        sys.exit(1)

    # Validate input files exist
    if not Path(chgcar_neutral).exists():
        print(f"Error: File not found: {chgcar_neutral}", file=sys.stderr)
        sys.exit(1)
    if not Path(chgcar_fukui).exists():
        print(f"Error: File not found: {chgcar_fukui}", file=sys.stderr)
        sys.exit(1)

    # Add FukuiGrid to Python path
    # Path: scripts/ -> fukui/ -> core/ -> teros/ -> external/FukuiGrid/
    script_dir = Path(__file__).parent
    fukui_grid_path = script_dir.parent.parent.parent / 'external' / 'FukuiGrid'

    if not fukui_grid_path.exists():
        print(f"Error: FukuiGrid not found at {fukui_grid_path}", file=sys.stderr)
        print("Please clone FukuiGrid to teros/external/FukuiGrid", file=sys.stderr)
        sys.exit(1)

    sys.path.insert(0, str(fukui_grid_path))

    # Import and run FukuiGrid
    try:
        import FukuiGrid
    except ImportError as e:
        print(f"Error importing FukuiGrid: {e}", file=sys.stderr)
        sys.exit(1)

    print("Running fukui_electrodes with:")
    print("  Neutral CHGCAR: {}".format(chgcar_neutral))
    print("  Fukui CHGCAR: {}".format(chgcar_fukui))
    print("  Epsilon: {}".format(epsilon))

    # Call FukuiGrid function
    # fukui_electrodes(FILE0, FILE1, Epsilon)
    # FILE0 = neutral charge density
    # FILE1 = Fukui function
    # Returns path to LOCPOT_FUKUI.vasp
    result = FukuiGrid.fukui_electrodes(chgcar_neutral, chgcar_fukui, epsilon)

    print(f"Output written to: {result}")


if __name__ == '__main__':
    main()
