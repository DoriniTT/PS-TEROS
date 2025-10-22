#!/usr/bin/env python
"""
Verification script for coordination_finder module.
Demonstrates finding coordinated O atoms for the surface Mn.
"""

from pathlib import Path
import numpy as np
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens


def main():
    print("=" * 70)
    print("Coordination Finder Verification")
    print("=" * 70)
    print()

    # Load structure
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    print(f"Loading structure from: {cif_path.name}")
    structure = load_structure(cif_path)
    print(f"Total atoms: {len(structure)}")
    print()

    # Find surface Mn
    print("Finding surface Mn atom...")
    mn_index, mn_position = find_surface_mn(structure)
    print(f"Surface Mn:")
    print(f"  Index: {mn_index}")
    print(f"  Position: [{mn_position[0]:.4f}, {mn_position[1]:.4f}, {mn_position[2]:.4f}]")
    print(f"  Z-coordinate: {mn_position[2]:.4f} Å")
    print()

    # Find coordinated O atoms
    print("Finding coordinated O atoms (cutoff = 2.5 Å)...")
    o_indices = find_coordinated_oxygens(structure, mn_index, cutoff=2.5)
    print(f"Found {len(o_indices)} coordinated O atoms:")
    print()

    # Display details for each coordinated O
    for i, o_idx in enumerate(o_indices, 1):
        o_position = structure[o_idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        print(f"  O#{i}:")
        print(f"    Index: {o_idx}")
        print(f"    Position: [{o_position[0]:.4f}, {o_position[1]:.4f}, {o_position[2]:.4f}]")
        print(f"    Distance from Mn: {distance:.4f} Å")
        print(f"    Z-coordinate: {o_position[2]:.4f} Å")
        print()

    # Summary
    print("=" * 70)
    print("Summary:")
    print(f"  - Surface Mn at index {mn_index} (z = {mn_position[2]:.4f} Å)")
    print(f"  - {len(o_indices)} coordinated O atoms within 2.5 Å")

    distances = [np.linalg.norm(mn_position - structure[idx].coords)
                 for idx in o_indices]
    print(f"  - Distance range: {min(distances):.4f} - {max(distances):.4f} Å")
    print(f"  - O atoms are sorted by distance from Mn")
    print("=" * 70)


if __name__ == "__main__":
    main()
