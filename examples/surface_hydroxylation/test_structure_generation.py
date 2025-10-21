#!/home/thiagotd/envs/aiida/bin/python
"""
Test script for surface_hydroxylation module - Structure Generation Only

This script tests ONLY the structure generation component of the
surface_hydroxylation module. It does NOT submit VASP jobs (too slow for testing).

Tests:
- Structure generation using generate_structures CalcFunction
- Hydroxylation mode (adding OH groups to O sites)
- Coverage-based deduplication
- Structure variant creation

Usage:
    source ~/envs/aiida/bin/activate
    python test_structure_generation.py
"""

import sys
from aiida import orm, load_profile
from ase.build import fcc111
from teros.core.surface_hydroxylation.tasks import generate_structures


def main():
    """Test structure generation component."""

    print("\n" + "="*70)
    print("SURFACE HYDROXYLATION - STRUCTURE GENERATION TEST")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded: psteros")

    # Create test slab (2x2 Pt(111) with O adlayer)
    print("\n2. Creating test structure...")
    print("   Building 2x2 Pt(111) slab with O adlayer")

    # Build Pt slab
    slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)
    print(f"   - Pt slab: {len(slab)} atoms")

    # Add O atoms on top of Pt surface (manually positioned above Pt atoms)
    # Get top layer Pt positions
    import numpy as np
    z_max = max(atom.position[2] for atom in slab)
    top_pt_atoms = [atom for atom in slab if abs(atom.position[2] - z_max) < 0.1]

    # Add O atoms 2.0 Å above each top Pt atom
    from ase import Atom
    for pt_atom in top_pt_atoms:
        o_pos = pt_atom.position.copy()
        o_pos[2] += 2.0  # Place O 2 Å above Pt
        slab.append(Atom('O', position=o_pos))

    slab.center(vacuum=10.0, axis=2)

    print(f"   - Total atoms: {len(slab)}")
    print(f"   - O atoms: {sum(1 for atom in slab if atom.symbol == 'O')}")

    # Convert to AiiDA
    structure = orm.StructureData(ase=slab)
    print("   ✓ Structure created")

    # Define surface generation parameters
    print("\n3. Setting up surface generation parameters...")
    surface_params = orm.Dict(dict={
        'mode': 'hydrogen',
        'species': 'O',
        'z_window': 0.5,
        'which_surface': 'top',
        'oh_dist': 0.98,
        'include_empty': False,
        'deduplicate_by_coverage': True,
        'coverage_bins': 3
    })

    print("   Parameters:")
    print(f"   - Mode: hydrogen (hydroxylation)")
    print(f"   - Species: O")
    print(f"   - Coverage bins: 3")
    print(f"   - Deduplication: enabled")

    # Run structure generation
    print("\n4. Running structure generation...")
    print("   This will create variants with different OH coverages")
    print("   (This is a CalcFunction - will be stored in provenance)")

    from aiida.engine import run
    result = run(generate_structures, structure=structure, params=surface_params)

    # Analyze results
    print("\n5. Analyzing results...")

    manifest = result['manifest'].get_dict()

    print(f"\n   ✓ Manifest created")
    print(f"   - Total variants generated: {len(manifest['variants'])}")

    # Show variant details
    print("\n   Variant details:")
    for idx, variant in enumerate(manifest['variants']):
        structure_key = f'structure_{idx}'
        if structure_key in result:
            variant_struct = result[structure_key]
            n_atoms = len(variant_struct.sites)
            coverage = variant.get('OH_coverage', 0.0)
            name = variant.get('name', 'unknown')

            print(f"   [{idx}] {name}")
            print(f"       - Coverage: {coverage:.2f}")
            print(f"       - Atoms: {n_atoms}")
            print(f"       - PK: {variant_struct.pk}")

    # Verify structures were created
    print("\n6. Verification...")

    expected_structures = len(manifest['variants'])
    actual_structures = sum(1 for key in result.keys() if key.startswith('structure_'))

    if expected_structures == actual_structures:
        print(f"   ✓ All {expected_structures} structures created successfully")
    else:
        print(f"   ✗ Mismatch: expected {expected_structures}, got {actual_structures}")
        return 1

    # Check coverage distribution
    coverages = [v.get('OH_coverage', 0.0) for v in manifest['variants']]
    print(f"   ✓ Coverage range: {min(coverages):.2f} - {max(coverages):.2f}")

    print("\n" + "="*70)
    print("STRUCTURE GENERATION TEST COMPLETED SUCCESSFULLY")
    print("="*70)
    print("\nNext steps:")
    print("  - Check provenance: verdi node show", result['manifest'].pk)
    print("  - View structure: verdi data core.structure show <PK>")
    print("\nTo test full workflow (with VASP):")
    print("  - Use test_full_workflow.py (requires VASP setup)")
    print("="*70 + "\n")

    return 0


if __name__ == '__main__':
    sys.exit(main())
