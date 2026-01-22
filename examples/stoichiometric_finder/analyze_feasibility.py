#!/usr/bin/env python
"""
Stoichiometric+Symmetric Surface Feasibility Analysis

This script demonstrates how to analyze which Miller indices have
stoichiometric AND symmetric surfaces before running expensive DFT calculations.

Usage:
    python analyze_feasibility.py

This is a pure Python analysis (no AiiDA needed) that helps decide:
- Which orientations can use the simple surface energy formula
- Which orientations need thermodynamics with chemical potentials
"""

import os
from pymatgen.core import Structure
from teros.core.surface_energy import (
    analyze_miller_feasibility,
    get_feasibility_summary,
    find_stoichiometric_symmetric_slabs,
    NoStoichiometricSymmetricSurfaceError,
)


def main():
    """Analyze feasibility for various materials."""

    print("\n" + "=" * 70)
    print("STOICHIOMETRIC+SYMMETRIC SURFACE FEASIBILITY ANALYSIS")
    print("=" * 70)

    # Example 1: FCC Gold (elemental metal - should be easy)
    print("\n" + "-" * 70)
    print("Example 1: FCC Gold (Au)")
    print("-" * 70)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    au_cif = os.path.join(script_dir, '../surface_thermo_gold/au.cif')

    if os.path.exists(au_cif):
        au_structure = Structure.from_file(au_cif)

        # Analyze common FCC orientations
        miller_indices = [(1, 1, 1), (1, 0, 0), (1, 1, 0)]

        reports = analyze_miller_feasibility(
            au_structure,
            miller_indices=miller_indices,
            min_slab_thickness=15.0,
        )

        print(get_feasibility_summary(reports))
    else:
        print(f"  [SKIPPED] Au CIF not found at: {au_cif}")

    # Example 2: B2 Intermetallic (PdIn) - more challenging
    print("\n" + "-" * 70)
    print("Example 2: B2 Intermetallic (PdIn)")
    print("-" * 70)

    # Create a simple B2 (CsCl-type) PdIn structure
    # Pd at (0,0,0), In at (0.5,0.5,0.5)
    from pymatgen.core import Lattice

    lattice = Lattice.cubic(3.25)  # Approximate lattice constant
    pdin_structure = Structure(
        lattice,
        ["Pd", "In"],
        [[0, 0, 0], [0.5, 0.5, 0.5]]
    )

    # B2 structures can have stoichiometric issues on certain surfaces
    miller_indices = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]

    reports = analyze_miller_feasibility(
        pdin_structure,
        miller_indices=miller_indices,
        min_slab_thickness=15.0,
    )

    print(get_feasibility_summary(reports))

    # Detailed analysis of a specific orientation
    print("\n  Detailed analysis for (1, 1, 0):")
    try:
        results = find_stoichiometric_symmetric_slabs(
            pdin_structure,
            miller_index=(1, 1, 0),
            min_slab_thickness=15.0,
            raise_if_not_found=True,
        )
        print(f"    Found {len(results)} valid slab(s):")
        for i, result in enumerate(results):
            print(f"      {i+1}. Thickness: {result.thickness_angstrom:.1f} A, "
                  f"Strategy: {result.strategy_used}")
    except NoStoichiometricSymmetricSurfaceError as e:
        print(f"    No valid surfaces found!")
        print(f"    Recommendation: Use thermodynamics module")

    # Example 3: Using thickness scan for difficult cases
    print("\n" + "-" * 70)
    print("Example 3: Thickness Scan Strategy")
    print("-" * 70)

    print("\n  Scanning thicknesses 10-30 A for PdIn (1,0,0):")

    results = find_stoichiometric_symmetric_slabs(
        pdin_structure,
        miller_index=(1, 0, 0),
        strategies=['thickness_scan'],  # Only use thickness scan
        raise_if_not_found=False,
    )

    if results:
        print(f"    Found {len(results)} valid slab(s) at different thicknesses:")
        for result in results[:5]:  # Show first 5
            print(f"      - {result.thickness_angstrom:.1f} A")
    else:
        print("    No valid surfaces found even with thickness scan")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("""
This analysis helps you decide:

1. For orientations with has_valid_surfaces=True:
   -> Use build_metal_surface_energy_workgraph() with simple formula
   -> Surface energy = (E_slab - N*E_bulk/atom) / 2A

2. For orientations with has_valid_surfaces=False:
   -> Use build_thermodynamics_workgraph() from teros.core.thermodynamics
   -> Surface energy = gamma(delta_mu) requires chemical potential

3. For mixed cases:
   -> Set require_stoichiometric_symmetric=True in workgraph
   -> Only valid terminations will be calculated
   -> Failed orientations will raise NoStoichiometricSymmetricSurfaceError
""")


if __name__ == '__main__':
    main()
