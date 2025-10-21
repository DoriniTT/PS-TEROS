# test_surface_detector.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn


def test_find_surface_mn_returns_topmost():
    """Test finding Mn with highest z-coordinate"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    mn_index, mn_position = find_surface_mn(structure)

    assert mn_index is not None
    assert mn_index >= 0
    assert len(mn_position) == 3  # x, y, z coordinates

    # Verify it's actually a Mn atom
    assert str(structure[mn_index].specie) == "Mn"

    # Verify it's the highest z among all Mn atoms
    mn_sites = [i for i, site in enumerate(structure) if str(site.specie) == "Mn"]
    for idx in mn_sites:
        assert structure[mn_index].coords[2] >= structure[idx].coords[2]


def test_find_surface_mn_returns_cartesian_coords():
    """Test that returned position is in Cartesian coordinates"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    mn_index, mn_position = find_surface_mn(structure)

    # Cartesian coordinates should match structure's cartesian coords
    import numpy as np
    expected_coords = structure[mn_index].coords
    assert np.allclose(mn_position, expected_coords, atol=1e-6)
