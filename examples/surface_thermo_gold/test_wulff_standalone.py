#!/usr/bin/env python
"""
Standalone Test for Wulff Shape Generation

This script tests the Wulff shape generation with mock surface energy data,
without needing to run actual VASP calculations. Use this to verify that:
1. The pymatgen WulffShape integration works
2. PDF, TXT, and PY outputs are generated correctly
3. The lowest-energy termination selection logic works

Usage:
    source ~/envs/aiida/bin/activate
    python test_wulff_standalone.py

This will create output files in the current directory:
- au_wulff_shape.pdf
- au_wulff_report.txt
- au_wulff_plot.py
"""

import os
import sys
import tempfile

# Test without AiiDA first to check pymatgen integration
def test_wulff_without_aiida():
    """Test Wulff shape generation using pymatgen directly."""
    print("\n" + "="*70)
    print("WULFF SHAPE STANDALONE TEST (without AiiDA)")
    print("="*70)

    try:
        from pymatgen.core.lattice import Lattice
        from pymatgen.analysis.wulff import WulffShape
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        print("✓ pymatgen and matplotlib imported successfully")
    except ImportError as e:
        print(f"✗ Import error: {e}")
        print("  Make sure pymatgen and matplotlib are installed")
        return False

    # Mock data for Au (FCC, a = 4.08 Å)
    print("\n1. Creating mock Au lattice (FCC, a = 4.08 Å)...")
    lattice = Lattice.cubic(4.08)
    print(f"   Lattice: a = {lattice.a:.3f} Å")

    # Mock surface energies (typical values for Au)
    miller_indices = [(1, 1, 1), (1, 0, 0), (1, 1, 0)]
    surface_energies = [0.79, 0.94, 1.33]  # J/m²

    print("\n2. Mock surface energies:")
    for miller, gamma in zip(miller_indices, surface_energies):
        print(f"   {miller}: {gamma:.2f} J/m²")

    # Create Wulff shape
    print("\n3. Creating WulffShape...")
    try:
        wulff = WulffShape(lattice, miller_indices, surface_energies)
        print("   ✓ WulffShape created successfully")
    except Exception as e:
        print(f"   ✗ Error creating WulffShape: {e}")
        return False

    # Print Wulff shape properties
    print("\n4. Wulff shape properties:")
    print(f"   Weighted surface energy: {wulff.weighted_surface_energy:.4f} J/m²")
    print(f"   Anisotropy: {wulff.anisotropy:.4f}")
    print(f"   Shape factor: {wulff.shape_factor:.4f}")
    print(f"   Volume: {wulff.volume:.6f}")
    print(f"   Effective radius: {wulff.effective_radius:.4f}")

    # Print area fractions
    print("\n5. Area fractions:")
    for miller, frac in wulff.area_fraction_dict.items():
        print(f"   {tuple(miller)}: {frac*100:.1f}%")

    # Generate plot
    print("\n6. Generating PDF plot...")
    try:
        fig = wulff.get_plot(
            color_set='RdBu',
            grid_off=True,
            axis_off=True,
            show_area=True,
            alpha=0.9,
        )
        fig.suptitle('Wulff Shape: Au (Test)', fontsize=14, fontweight='bold')

        output_pdf = 'au_wulff_shape_test.pdf'
        fig.savefig(output_pdf, format='pdf', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"   ✓ Saved: {output_pdf}")
    except Exception as e:
        print(f"   ✗ Error generating plot: {e}")
        return False

    print("\n" + "="*70)
    print("✓ STANDALONE TEST PASSED")
    print("="*70)
    return True


def test_wulff_with_aiida():
    """Test the full AiiDA calcfunction with mock data."""
    print("\n" + "="*70)
    print("WULFF SHAPE TEST WITH AIIDA CALCFUNCTION")
    print("="*70)

    try:
        from aiida import load_profile, orm
        load_profile('presto')
        print("✓ AiiDA profile loaded")
    except Exception as e:
        print(f"✗ Could not load AiiDA: {e}")
        return False

    try:
        from teros.core.surface_energy.wulff import (
            generate_wulff_shape_data,
            select_lowest_energy_terminations,
        )
        print("✓ Wulff module imported")
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False

    # Create mock surface_energies Dict (simulating workflow output)
    print("\n1. Creating mock surface_energies Dict...")
    mock_surface_energies = {
        'hkl_111': {
            'term_0': {
                'gamma_J_m2': 0.79,
                'gamma_eV_A2': 0.0493,
                'area_A2': 21.5,
                'E_slab_eV': -145.23,
                'E_bulk_per_atom_eV': -3.81,
                'N_slab': 36,
                'N_bulk': 4,
                'compound_type': 'elemental',
                'formula': 'Au',
                'elements': ['Au'],
                'composition': {'Au': 4},
                'slab_composition': {'Au': 36},
                'is_stoichiometric': True,
            },
            'term_1': {
                'gamma_J_m2': 0.85,  # Higher energy termination
                'gamma_eV_A2': 0.0531,
                'area_A2': 21.5,
                'E_slab_eV': -144.95,
                'E_bulk_per_atom_eV': -3.81,
                'N_slab': 36,
                'N_bulk': 4,
                'compound_type': 'elemental',
                'formula': 'Au',
                'elements': ['Au'],
                'composition': {'Au': 4},
                'slab_composition': {'Au': 36},
                'is_stoichiometric': True,
            },
        },
        'hkl_100': {
            'term_0': {
                'gamma_J_m2': 0.94,
                'gamma_eV_A2': 0.0587,
                'area_A2': 16.6,
                'E_slab_eV': -140.12,
                'E_bulk_per_atom_eV': -3.81,
                'N_slab': 36,
                'N_bulk': 4,
                'compound_type': 'elemental',
                'formula': 'Au',
                'elements': ['Au'],
                'composition': {'Au': 4},
                'slab_composition': {'Au': 36},
                'is_stoichiometric': True,
            },
        },
        'hkl_110': {
            'term_0': {
                'gamma_J_m2': 1.33,
                'gamma_eV_A2': 0.0831,
                'area_A2': 23.5,
                'E_slab_eV': -138.56,
                'E_bulk_per_atom_eV': -3.81,
                'N_slab': 36,
                'N_bulk': 4,
                'compound_type': 'elemental',
                'formula': 'Au',
                'elements': ['Au'],
                'composition': {'Au': 4},
                'slab_composition': {'Au': 36},
                'is_stoichiometric': True,
            },
        },
    }

    # Test select_lowest_energy_terminations
    print("\n2. Testing termination selection...")
    selected = select_lowest_energy_terminations(mock_surface_energies)
    print(f"   Selected terminations:")
    for miller, data in selected.items():
        print(f"     {miller}: {data['termination']} (γ = {data['gamma_J_m2']:.2f} J/m²)")

    # Verify (111) picked the lower energy term_0
    assert selected[(1, 1, 1)]['termination'] == 'term_0', "Should select term_0 for (111)"
    assert selected[(1, 1, 1)]['gamma_J_m2'] == 0.79, "Wrong energy for (111)"
    print("   ✓ Termination selection correct")

    # Create mock bulk structure (FCC Au)
    print("\n3. Creating mock bulk structure...")
    from ase.build import bulk
    au_bulk = bulk('Au', 'fcc', a=4.08)
    bulk_structure = orm.StructureData(ase=au_bulk)
    print(f"   Structure: {bulk_structure.get_formula()}")
    print(f"   Cell: {au_bulk.cell[0][0]:.3f} Å")

    # Create AiiDA Dict node
    surface_energies_node = orm.Dict(dict=mock_surface_energies)

    # Run the calcfunction
    print("\n4. Running generate_wulff_shape_data calcfunction...")
    try:
        result = generate_wulff_shape_data(
            surface_energies=surface_energies_node,
            bulk_structure=bulk_structure,
        )
        print("   ✓ Calcfunction completed")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

    # Check outputs
    print("\n5. Checking outputs...")

    # Check wulff_data
    wulff_data = result['wulff_data'].get_dict()
    print(f"   wulff_data keys: {list(wulff_data.keys())}")
    print(f"   Weighted surface energy: {wulff_data['weighted_surface_energy_J_m2']:.4f} J/m²")
    print(f"   Shape factor: {wulff_data['shape_factor']:.4f}")
    print(f"   Anisotropy: {wulff_data['anisotropy']:.4f}")

    # Check file outputs
    print(f"\n   wulff_plot_pdf: PK={result['wulff_plot_pdf'].pk}")
    print(f"   wulff_report_txt: PK={result['wulff_report_txt'].pk}")
    print(f"   wulff_script_py: PK={result['wulff_script_py'].pk}")

    # Extract files to current directory
    print("\n6. Extracting output files...")
    import shutil

    for name, node in [
        ('au_wulff_shape.pdf', result['wulff_plot_pdf']),
        ('au_wulff_report.txt', result['wulff_report_txt']),
        ('au_wulff_plot.py', result['wulff_script_py']),
    ]:
        with node.open(mode='rb') as f:
            with open(name, 'wb') as out:
                out.write(f.read())
        print(f"   ✓ Saved: {name}")

    print("\n" + "="*70)
    print("✓ AIIDA CALCFUNCTION TEST PASSED")
    print("="*70)

    # Show report contents
    print("\n--- WULFF REPORT PREVIEW ---")
    with open('au_wulff_report.txt', 'r') as f:
        print(f.read())

    return True


if __name__ == '__main__':
    print("Testing Wulff Shape Generation")
    print("="*70)

    # First test without AiiDA (just pymatgen)
    standalone_ok = test_wulff_without_aiida()

    if not standalone_ok:
        print("\n✗ Standalone test failed. Fix pymatgen issues first.")
        sys.exit(1)

    # Then test with AiiDA
    print("\n" + "="*70)
    print("Now testing with AiiDA calcfunction...")
    print("="*70)

    aiida_ok = test_wulff_with_aiida()

    if aiida_ok:
        print("\n" + "="*70)
        print("ALL TESTS PASSED!")
        print("="*70)
        print("\nOutput files created:")
        print("  - au_wulff_shape_test.pdf (from standalone test)")
        print("  - au_wulff_shape.pdf (from AiiDA test)")
        print("  - au_wulff_report.txt")
        print("  - au_wulff_plot.py")
        print("\nYou can now run the standalone script to customize the plot:")
        print("  python au_wulff_plot.py")
    else:
        print("\n✗ AiiDA test failed")
        sys.exit(1)
