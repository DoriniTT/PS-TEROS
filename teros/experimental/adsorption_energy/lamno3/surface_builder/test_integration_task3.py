#!/usr/bin/env python
"""Integration test for Task 3: Find coordinated O atoms"""

from pathlib import Path
import numpy as np
from structure_loader import load_structure, validate_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens


def main():
    print("=" * 70)
    print("Task 3 Integration Test: Find Coordinated O Atoms")
    print("=" * 70)
    print()

    # Load and validate structure
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    validation = validate_structure(structure)

    assert validation["valid"], "Structure should be valid"
    assert validation["num_mn"] > 0, "Should have Mn atoms"
    assert validation["num_o"] > 0, "Should have O atoms"
    print("✓ Structure loaded and validated")

    # Find surface Mn
    mn_index, mn_position = find_surface_mn(structure)
    assert mn_index >= 0, "Should find surface Mn"
    assert str(structure[mn_index].specie) == "Mn", "Should be Mn atom"
    print(f"✓ Surface Mn found at index {mn_index}")

    # Find coordinated O atoms
    o_indices = find_coordinated_oxygens(structure, mn_index, cutoff=2.5)
    assert isinstance(o_indices, list), "Should return list"
    assert len(o_indices) > 0, "Should find coordinated O atoms"
    print(f"✓ Found {len(o_indices)} coordinated O atoms")

    # Verify all are O atoms
    for idx in o_indices:
        assert str(structure[idx].specie) == "O", f"Index {idx} should be O atom"
    print("✓ All coordinated atoms are O")

    # Verify distances
    distances = []
    for idx in o_indices:
        o_position = structure[idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        distances.append(distance)
        assert distance <= 2.5, f"Distance {distance:.3f} exceeds cutoff"
        assert distance > 1.5, f"Distance {distance:.3f} unreasonably small"
    print(f"✓ All distances within range: {min(distances):.3f}-{max(distances):.3f} Å")

    # Verify sorting
    assert distances == sorted(distances), "O atoms should be sorted by distance"
    print("✓ O atoms sorted by distance from Mn")

    print()
    print("=" * 70)
    print("Task 3 Integration Test: PASSED")
    print("=" * 70)
    print()
    print("Summary:")
    print(f"  - Loaded structure with {len(structure)} atoms")
    print(f"  - Found surface Mn at index {mn_index}")
    print(f"  - Found {len(o_indices)} coordinated O atoms")
    print(f"  - Coordinated O indices: {o_indices}")
    print(f"  - Mn-O distances: {[f'{d:.3f} Å' for d in distances]}")
    print()


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print()
        print("=" * 70)
        print(f"FAILED: {e}")
        print("=" * 70)
        exit(1)
    except Exception as e:
        print()
        print("=" * 70)
        print(f"ERROR: {e}")
        print("=" * 70)
        exit(1)
