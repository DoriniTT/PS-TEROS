# test_structure_loader.py
import pytest
from pathlib import Path
from structure_loader import load_structure, validate_structure


def test_load_structure_from_cif():
    """Test loading CIF file returns valid Structure object"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    assert structure is not None
    assert len(structure) > 0  # Has atoms
    assert structure.lattice is not None


def test_validate_structure_has_required_elements():
    """Test validation checks for Mn and O elements"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)

    result = validate_structure(structure)

    assert result["valid"] is True
    assert "Mn" in result["elements"]
    assert "O" in result["elements"]
    assert result["num_mn"] > 0
    assert result["num_o"] > 0


def test_load_structure_invalid_path():
    """Test loading non-existent file raises error"""
    with pytest.raises(FileNotFoundError):
        load_structure("nonexistent.cif")
