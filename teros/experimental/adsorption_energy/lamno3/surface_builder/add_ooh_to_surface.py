#!/usr/bin/env python
"""
Add OOH radical to LaMnO3 surface structure.

This script:
1. Loads a surface slab structure (CIF)
2. Creates a 2x2x1 supercell
3. Identifies the topmost Mn atom (surface Mn)
4. Finds O atoms coordinated to this Mn
5. Selects the most exposed O (highest z)
6. Adds O and H to form OOH radical
7. Exports modified structure to new CIF file
"""

import sys
from pathlib import Path
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

from structure_loader import load_structure, validate_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen
from ooh_constructor import construct_ooh_radical
from structure_exporter import export_structure_to_cif


def main():
    """Main execution function"""

    # Configuration
    input_cif = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    output_cif = Path(__file__).parent / "LaMnO3_100_A4_surface_2x2x1_OOH.cif"

    print("=" * 60)
    print("OOH Surface Modification Script")
    print("=" * 60)
    print()

    # Step 1: Load and validate structure
    print(f"Loading structure from: {input_cif}")
    structure = load_structure(input_cif)

    validation = validate_structure(structure)
    if not validation["valid"]:
        print("ERROR: Structure validation failed")
        print(f"  Missing required elements")
        sys.exit(1)

    print(f"✓ Structure loaded successfully")
    print(f"  Total atoms: {len(structure)}")
    print(f"  Elements: {', '.join(validation['elements'])}")
    print(f"  Mn atoms: {validation['num_mn']}")
    print(f"  O atoms: {validation['num_o']}")
    print()

    # Step 2: Create 2x2x1 supercell
    print("Creating 2x2x1 supercell...")
    structure.make_supercell([2, 2, 1])
    print(f"✓ Supercell created")
    print(f"  New total atoms: {len(structure)}")
    print(f"  New Mn atoms: {sum(1 for site in structure if str(site.specie) == 'Mn')}")
    print(f"  New O atoms: {sum(1 for site in structure if str(site.specie) == 'O')}")
    print()

    # Step 3: Find surface Mn
    print("Finding surface Mn atom...")
    mn_index, mn_position = find_surface_mn(structure)
    print(f"✓ Surface Mn identified")
    print(f"  Index: {mn_index}")
    print(f"  Position (Cartesian): [{mn_position[0]:.4f}, {mn_position[1]:.4f}, {mn_position[2]:.4f}]")
    print()

    # Step 4: Find coordinated O atoms
    print("Finding coordinated O atoms...")
    o_indices = find_coordinated_oxygens(structure, mn_index, cutoff=2.5)
    print(f"✓ Found {len(o_indices)} coordinated O atoms")
    for i, o_idx in enumerate(o_indices, 1):
        o_pos = structure[o_idx].coords
        distance = np.linalg.norm(mn_position - o_pos)
        print(f"  O{i}: index={o_idx}, distance={distance:.3f} Å")
    print()

    # Step 5: Select most exposed O
    print("Selecting most exposed O atom...")
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)
    print(f"✓ Most exposed O selected")
    print(f"  Index: {exposed_o_idx}")
    print(f"  Position (Cartesian): [{exposed_o_pos[0]:.4f}, {exposed_o_pos[1]:.4f}, {exposed_o_pos[2]:.4f}]")
    print(f"  Z-coordinate: {exposed_o_pos[2]:.4f} Å")
    print()

    # Step 6: Convert to ASE and construct OOH
    print("Constructing OOH radical...")
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    print(f"✓ OOH radical constructed")
    print(f"  Original atoms: {len(ase_atoms)}")
    print(f"  Modified atoms: {len(modified_atoms)}")
    print(f"  Added: 1 O + 1 H")

    # Get positions of new atoms
    positions = modified_atoms.get_positions()
    new_o_pos = positions[len(ase_atoms)]
    h_pos = positions[len(ase_atoms) + 1]

    oo_distance = np.linalg.norm(new_o_pos - exposed_o_pos)
    oh_distance = np.linalg.norm(h_pos - new_o_pos)

    print(f"  O-O bond length: {oo_distance:.3f} Å")
    print(f"  O-H bond length: {oh_distance:.3f} Å")
    print()

    # Step 7: Export modified structure
    print(f"Exporting modified structure to: {output_cif}")
    export_structure_to_cif(modified_atoms, output_cif)
    print(f"✓ Structure exported successfully")
    print()

    print("=" * 60)
    print("OOH modification completed successfully!")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
