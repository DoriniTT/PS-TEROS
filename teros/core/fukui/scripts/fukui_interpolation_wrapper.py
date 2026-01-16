#!/usr/bin/env python
"""Wrapper script for FukuiGrid interpolation (used by aiida-shell).

This script is called by aiida-shell with CHGCAR files in the current directory.
It imports FukuiGrid using an absolute path and runs the interpolation.

Usage:
    python fukui_interpolation_wrapper.py <fukui_type> <delta_n_csv> <file1> <file2> <file3> <file4>

Example:
    python fukui_interpolation_wrapper.py plus 0.0,0.05,0.10,0.15 CHGCAR_0.00 CHGCAR_0.05 CHGCAR_0.10 CHGCAR_0.15

Output:
    CHGCAR_FUKUI.vasp - Fukui function mapped onto charge density grid
"""

import sys
import os
import numpy as np
from pathlib import Path

# Compute FukuiGrid path relative to this script's location
# Script is at: teros/core/fukui/scripts/fukui_interpolation_wrapper.py
# FukuiGrid is at: teros/external/FukuiGrid/
_SCRIPT_DIR = Path(__file__).resolve().parent
FUKUI_GRID_PATH = str(_SCRIPT_DIR.parent.parent.parent / 'external' / 'FukuiGrid')


def main():
    """Run FukuiGrid interpolation with command-line arguments."""
    # Parse arguments: fukui_type, delta_n values (comma-separated), file1, file2, file3, file4
    if len(sys.argv) < 7:
        print("Usage: fukui_interpolation_wrapper.py <fukui_type> <delta_n_csv> <file1> <file2> <file3> <file4>")
        print("  fukui_type: 'plus' (nucleophilic f+) or 'minus' (electrophilic f-)")
        print("  delta_n_csv: comma-separated delta_n values (e.g., '0.0,0.05,0.10,0.15')")
        print("  file1-4: CHGCAR files corresponding to each delta_n value")
        sys.exit(1)

    fukui_type = sys.argv[1]  # 'plus' or 'minus'
    delta_n_csv = sys.argv[2]  # e.g., "0.0,0.05,0.10,0.15"
    files = sys.argv[3:7]

    # Validate fukui_type
    if fukui_type not in ('plus', 'minus'):
        print(f"ERROR: fukui_type must be 'plus' or 'minus', got '{fukui_type}'")
        sys.exit(1)

    # Parse delta_n values
    try:
        delta_n = np.array([float(x) for x in delta_n_csv.split(',')])
    except ValueError as e:
        print(f"ERROR: Could not parse delta_n values: {e}")
        sys.exit(1)

    if len(delta_n) != len(files):
        print(f"ERROR: Number of delta_n values ({len(delta_n)}) must match number of files ({len(files)})")
        sys.exit(1)

    # Verify all files exist
    for f in files:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            sys.exit(1)

    print(f"FukuiGrid Interpolation Wrapper")
    print(f"  Fukui type: {fukui_type}")
    print(f"  Delta N values: {delta_n.tolist()}")
    print(f"  Input files: {files}")

    # Sort by delta_n descending (required by FukuiGrid)
    sorted_indices = np.argsort(delta_n)[::-1]
    sorted_files = [files[i] for i in sorted_indices]
    sorted_dn = delta_n[sorted_indices]

    # Adjust sign for fukui_type
    # For f+: electrons are removed, so we use positive delta_n
    # For f-: electrons are added, so we use negative delta_n
    if fukui_type == 'minus':
        sorted_dn = -sorted_dn

    print(f"  Sorted files (by delta_n descending): {sorted_files}")
    print(f"  Sorted delta_n for FukuiGrid: {sorted_dn.tolist()}")

    # Import FukuiGrid from absolute path
    sys.path.insert(0, FUKUI_GRID_PATH)
    try:
        from FukuiGrid import Fukui_interpolation
    except ImportError as e:
        print(f"ERROR: Could not import FukuiGrid: {e}")
        print(f"  Expected location: {FUKUI_GRID_PATH}")
        sys.exit(1)

    # Run interpolation (files are in current working directory)
    print(f"\nRunning Fukui_interpolation...")
    Fukui_interpolation(sorted_files[0], sorted_files[1], sorted_files[2], sorted_files[3], dn=sorted_dn)

    # Verify output was created
    output_file = 'CHGCAR_FUKUI.vasp'
    if os.path.exists(output_file):
        size_mb = os.path.getsize(output_file) / (1024 * 1024)
        print(f"\nSUCCESS: {output_file} generated ({size_mb:.1f} MB)")
    else:
        print(f"\nERROR: Expected output file '{output_file}' was not created")
        sys.exit(1)


if __name__ == '__main__':
    main()
