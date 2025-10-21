# test_oxygen_selector.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen


def test_select_most_exposed_oxygen_returns_highest_z():
    """Test selecting O with highest z-coordinate from coordinated O atoms"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)

    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    assert exposed_o_idx in o_indices
    assert len(exposed_o_pos) == 3

    # Verify it's the highest z among coordinated O
    for idx in o_indices:
        assert exposed_o_pos[2] >= structure[idx].coords[2]


def test_select_most_exposed_oxygen_is_actually_oxygen():
    """Test that selected atom is actually O"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)

    exposed_o_idx, _ = select_most_exposed_oxygen(structure, o_indices)

    assert str(structure[exposed_o_idx].specie) == "O"


def test_select_most_exposed_oxygen_empty_list_raises_error():
    """Test that empty O list raises appropriate error"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    with pytest.raises(ValueError, match="No oxygen atoms provided"):
        select_most_exposed_oxygen(structure, [])
