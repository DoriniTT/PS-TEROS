#!/usr/bin/env python
"""Wrapper script for FukuiGrid.Perturbative_point() function.

This script is called by the run_perturbative_expansion_calcfunc to compute
the perturbative model potential using the c-DFT expansion:

    ΔU(r) = q·Φ(r) - q·ΔN·vf±(r)

Where:
    Φ(r) = electrostatic potential (LOCPOT from neutral slab)
    vf±(r) = Fukui potential (LOCPOT_FUKUI.vasp from Phase 2)
    q = charge of the probe (point charge model)
    ΔN = electron transfer

Usage:
    python perturbative_expansion_wrapper.py <locpot> <fukui_locpot> <q> <dN>

Arguments:
    locpot: Path to LOCPOT file (electrostatic potential of neutral slab)
    fukui_locpot: Path to Fukui potential file (e.g., LOCPOT_FUKUI.vasp)
    q: Charge of probe in |e|
    dN: Electron transfer ΔN

Output:
    MODELPOT_LOCPOT.vasp in the current working directory
"""

import sys
from pathlib import Path


def main():
    """Run FukuiGrid Perturbative_point with command-line arguments."""
    if len(sys.argv) != 5:
        print(
            "Usage: perturbative_expansion_wrapper.py <locpot> <fukui_locpot> <q> <dN>",
            file=sys.stderr,
        )
        print("\nArguments:", file=sys.stderr)
        print("  locpot: Path to LOCPOT file (neutral slab)", file=sys.stderr)
        print("  fukui_locpot: Path to Fukui potential file", file=sys.stderr)
        print("  q: Charge of probe in |e|", file=sys.stderr)
        print("  dN: Electron transfer", file=sys.stderr)
        sys.exit(1)

    locpot_file = sys.argv[1]
    fukui_locpot_file = sys.argv[2]

    try:
        q = float(sys.argv[3])
    except ValueError:
        print(f"Error: q must be a number, got '{sys.argv[3]}'", file=sys.stderr)
        sys.exit(1)

    try:
        dN = float(sys.argv[4])
    except ValueError:
        print(f"Error: dN must be a number, got '{sys.argv[4]}'", file=sys.stderr)
        sys.exit(1)

    # Validate input files exist
    if not Path(locpot_file).exists():
        print(f"Error: File not found: {locpot_file}", file=sys.stderr)
        sys.exit(1)
    if not Path(fukui_locpot_file).exists():
        print(f"Error: File not found: {fukui_locpot_file}", file=sys.stderr)
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

    print("Running Perturbative_point with:")
    print("  LOCPOT (neutral): {}".format(locpot_file))
    print("  Fukui potential: {}".format(fukui_locpot_file))
    print("  Probe charge q: {}".format(q))
    print("  Electron transfer dN: {}".format(dN))

    # Call FukuiGrid function
    # Perturbative_point(FILE1, FILE2, q, N)
    # FILE1 = electrostatic potential (LOCPOT)
    # FILE2 = Fukui potential (LOCPOT_FUKUI.vasp)
    # q = probe charge
    # N = electron transfer
    # Returns path to MODELPOT_LOCPOT.vasp
    result = FukuiGrid.Perturbative_point(locpot_file, fukui_locpot_file, q, dN)

    print(f"Output written to: {result}")


if __name__ == '__main__':
    main()
