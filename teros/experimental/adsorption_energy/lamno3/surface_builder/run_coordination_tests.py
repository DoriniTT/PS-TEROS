#!/usr/bin/env python
"""Simple test runner for coordination_finder without pytest"""

from pathlib import Path
import numpy as np
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens


def test_find_coordinated_oxygens_returns_list():
    """Test finding O atoms coordinated to surface Mn"""
    print("Test 1: find_coordinated_oxygens_returns_list...")
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index)

    assert isinstance(o_indices, list), "Result should be a list"
    assert len(o_indices) > 0, "Should have at least some coordinated O"

    # Verify all are actually O atoms
    for idx in o_indices:
        assert str(structure[idx].specie) == "O", f"Index {idx} is not an O atom"

    print(f"  PASS - Found {len(o_indices)} coordinated O atoms")


def test_coordinated_oxygens_within_distance():
    """Test that all coordinated O are within reasonable Mn-O bond distance"""
    print("Test 2: coordinated_oxygens_within_distance...")
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, mn_position = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index, cutoff=2.5)

    for idx in o_indices:
        o_position = structure[idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        assert distance <= 2.5, f"O atom {idx} is beyond cutoff: {distance:.3f} Å"
        assert distance > 1.5, f"O atom {idx} is too close: {distance:.3f} Å"

    print(f"  PASS - All {len(o_indices)} O atoms within distance cutoff")


def test_coordinated_oxygens_sorted_by_distance():
    """Test that returned O atoms are sorted by distance to Mn"""
    print("Test 3: coordinated_oxygens_sorted_by_distance...")
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, mn_position = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index)

    distances = []
    for idx in o_indices:
        o_position = structure[idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        distances.append(distance)

    # Check sorted
    assert distances == sorted(distances), "Distances are not sorted"

    print(f"  PASS - O atoms correctly sorted by distance")
    print(f"  Distances: {[f'{d:.3f}' for d in distances]}")


if __name__ == "__main__":
    print("=" * 60)
    print("Running Coordination Finder Tests")
    print("=" * 60)
    print()

    try:
        test_find_coordinated_oxygens_returns_list()
        test_coordinated_oxygens_within_distance()
        test_coordinated_oxygens_sorted_by_distance()

        print()
        print("=" * 60)
        print("All tests PASSED!")
        print("=" * 60)
    except AssertionError as e:
        print()
        print("=" * 60)
        print(f"Test FAILED: {e}")
        print("=" * 60)
        exit(1)
    except Exception as e:
        print()
        print("=" * 60)
        print(f"Test ERROR: {e}")
        print("=" * 60)
        exit(1)
