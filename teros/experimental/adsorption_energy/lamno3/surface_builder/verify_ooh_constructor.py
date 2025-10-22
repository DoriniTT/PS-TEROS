#!/usr/bin/env python
"""
Verify OOH constructor implementation with detailed output.
"""

from pathlib import Path
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen
from ooh_constructor import construct_ooh_radical


def main():
    print("=" * 70)
    print("OOH Constructor Verification")
    print("=" * 70)
    print()

    # Load structure
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    print(f"Loading structure from: {cif_path}")
    structure = load_structure(cif_path)
    print(f"Original structure has {len(structure)} atoms")
    print()

    # Find surface Mn
    print("Finding surface Mn atom...")
    mn_index, mn_position = find_surface_mn(structure)
    print(f"Surface Mn index: {mn_index}")
    print(f"Surface Mn position: [{mn_position[0]:.4f}, {mn_position[1]:.4f}, {mn_position[2]:.4f}]")
    print()

    # Find coordinated O atoms
    print("Finding coordinated O atoms...")
    o_indices = find_coordinated_oxygens(structure, mn_index)
    print(f"Found {len(o_indices)} coordinated O atoms")
    print()

    # Select most exposed O
    print("Selecting most exposed O atom...")
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)
    print(f"Most exposed O index: {exposed_o_idx}")
    print(f"Most exposed O position: [{exposed_o_pos[0]:.4f}, {exposed_o_pos[1]:.4f}, {exposed_o_pos[2]:.4f}]")
    print(f"Z-coordinate: {exposed_o_pos[2]:.4f} Å")
    print()

    # Convert to ASE
    print("Converting to ASE Atoms...")
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    print(f"ASE structure has {len(ase_atoms)} atoms")
    print()

    # Construct OOH radical
    print("Constructing OOH radical...")
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)
    print(f"Modified structure has {len(modified_atoms)} atoms")
    print(f"Added {len(modified_atoms) - len(ase_atoms)} atoms")
    print()

    # Get positions of new atoms
    positions = modified_atoms.get_positions()
    symbols = modified_atoms.get_chemical_symbols()

    new_o_idx = len(ase_atoms)
    h_idx = len(ase_atoms) + 1

    new_o_pos = positions[new_o_idx]
    h_pos = positions[h_idx]

    print("New atoms added:")
    print(f"  Atom {new_o_idx}: {symbols[new_o_idx]} at [{new_o_pos[0]:.4f}, {new_o_pos[1]:.4f}, {new_o_pos[2]:.4f}]")
    print(f"  Atom {h_idx}: {symbols[h_idx]} at [{h_pos[0]:.4f}, {h_pos[1]:.4f}, {h_pos[2]:.4f}]")
    print()

    # Calculate and verify bond lengths
    print("Verifying bond lengths...")
    oo_distance = np.linalg.norm(new_o_pos - exposed_o_pos)
    oh_distance = np.linalg.norm(h_pos - new_o_pos)

    print(f"  O-O bond length: {oo_distance:.4f} Å (expected: ~1.45 Å)")
    print(f"  O-H bond length: {oh_distance:.4f} Å (expected: ~0.97 Å)")
    print()

    # Verify orientation
    print("Verifying orientation...")
    print(f"  Base O z: {exposed_o_pos[2]:.4f} Å")
    print(f"  New O z: {new_o_pos[2]:.4f} Å (delta: +{new_o_pos[2] - exposed_o_pos[2]:.4f} Å)")
    print(f"  H z: {h_pos[2]:.4f} Å (delta from new O: +{h_pos[2] - new_o_pos[2]:.4f} Å)")
    print()

    # Verify tests
    print("Verification checks:")
    checks = [
        (len(modified_atoms) == len(ase_atoms) + 2, "Added exactly 2 atoms"),
        (symbols[new_o_idx] == 'O', "First new atom is O"),
        (symbols[h_idx] == 'H', "Second new atom is H"),
        (1.3 < oo_distance < 1.6, "O-O bond length in valid range (1.3-1.6 Å)"),
        (0.85 < oh_distance < 1.1, "O-H bond length in valid range (0.85-1.1 Å)"),
        (new_o_pos[2] > exposed_o_pos[2], "New O is above base O"),
        (h_pos[2] > new_o_pos[2], "H is above new O"),
    ]

    all_passed = True
    for passed, description in checks:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {description}")
        if not passed:
            all_passed = False

    print()
    print("=" * 70)
    if all_passed:
        print("SUCCESS: All verification checks passed!")
    else:
        print("FAILURE: Some verification checks failed!")
    print("=" * 70)

    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())
