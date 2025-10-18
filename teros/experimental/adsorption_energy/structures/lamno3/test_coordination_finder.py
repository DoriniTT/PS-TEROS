# test_coordination_finder.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens


def test_find_coordinated_oxygens_returns_list():
    """Test finding O atoms coordinated to surface Mn"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index)

    assert isinstance(o_indices, list)
    assert len(o_indices) > 0  # Should have at least some coordinated O

    # Verify all are actually O atoms
    for idx in o_indices:
        assert str(structure[idx].specie) == "O"


def test_coordinated_oxygens_within_distance():
    """Test that all coordinated O are within reasonable Mn-O bond distance"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, mn_position = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index, cutoff=2.5)

    import numpy as np
    for idx in o_indices:
        o_position = structure[idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        assert distance <= 2.5  # Within cutoff
        assert distance > 1.5  # Not unreasonably close


def test_coordinated_oxygens_sorted_by_distance():
    """Test that returned O atoms are sorted by distance to Mn"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, mn_position = find_surface_mn(structure)

    o_indices = find_coordinated_oxygens(structure, mn_index)

    import numpy as np
    distances = []
    for idx in o_indices:
        o_position = structure[idx].coords
        distance = np.linalg.norm(mn_position - o_position)
        distances.append(distance)

    # Check sorted
    assert distances == sorted(distances)
