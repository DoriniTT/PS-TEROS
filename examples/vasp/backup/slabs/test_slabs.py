#!/usr/bin/env python
"""
Simple test script to verify slab generation works correctly.
This tests just the slab generation function without running a full workflow.
"""

from aiida import load_profile, orm
from ase.io import read
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def test_slab_generation():
    """Test slab generation with a simple structure."""

    # Load AiiDA profile
    print("Loading AiiDA profile...")
    load_profile()

    # Load a test structure
    structures_dir = '/home/thiagotd/git/PS-TEROS/teros/structures'
    structure_file = f'{structures_dir}/ag3po4.cif'

    print(f"Loading structure from: {structure_file}")
    atoms = read(structure_file)

    print(f"Structure: {atoms.get_chemical_formula()}")
    print(f"Number of atoms: {len(atoms)}")

    # Test slab generation parameters
    miller_indices = [1, 0, 0]
    min_slab_thickness = 10.0
    min_vacuum_thickness = 15.0

    print(f"\nGenerating slabs for Miller indices: {miller_indices}")
    print(f"Min slab thickness: {min_slab_thickness} Å")
    print(f"Min vacuum thickness: {min_vacuum_thickness} Å")

    # Call slab generation directly using pymatgen
    print("\nGenerating slabs using Pymatgen...")

    # Convert to pymatgen structure
    adaptor = AseAtomsAdaptor()
    bulk_structure = adaptor.get_structure(atoms)

    # Get primitive cell if requested
    primitive = True
    if primitive:
        analyzer = SpacegroupAnalyzer(bulk_structure)
        bulk_structure = analyzer.get_primitive_standard_structure()
        print(f"Reduced to primitive cell: {bulk_structure.composition.reduced_formula}")

    # Create slab generator
    py_miller_indices = tuple(miller_indices)
    slab_gen = SlabGenerator(
        bulk_structure,
        py_miller_indices,
        min_slab_thickness,
        min_vacuum_thickness,
        lll_reduce=False,
        center_slab=True,
        max_normal_search=None,
        in_unit_planes=False
    )

    # Generate slabs
    slabs = slab_gen.get_slabs(symmetrize=False)
    print(f"\n✓ Generated {len(slabs)} slab terminations")

    # Convert to AiiDA structures
    result_slabs = {}
    for i, slab in enumerate(slabs):
        ortho_slab = slab.get_orthogonal_c_slab()
        super_slab = ortho_slab.make_supercell((1, 1, 1))
        ase_atoms = adaptor.get_atoms(super_slab)
        result_slabs[f"term_{i}"] = orm.StructureData(ase=ase_atoms)

    print(f"Termination identifiers: {list(result_slabs.keys())}")

    # Check each slab
    for term_id, slab_struct in result_slabs.items():
        slab_atoms = slab_struct.get_ase()
        print(f"\n{term_id}:")
        print(f"  Formula: {slab_atoms.get_chemical_formula()}")
        print(f"  Number of atoms: {len(slab_atoms)}")
        print(f"  Cell (Å):")
        for i, vec in enumerate(slab_atoms.get_cell()):
            print(f"    [{i}] {vec[0]:8.3f} {vec[1]:8.3f} {vec[2]:8.3f}")

    print("\n✓ All tests passed!")
    return {'slabs': result_slabs}

if __name__ == '__main__':
    try:
        result = test_slab_generation()
    except Exception as e:
        print(f"\n✗ Test failed!")
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        import sys
        sys.exit(1)
