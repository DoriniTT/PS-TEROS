# test_structure_exporter.py
import pytest
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen
from ooh_constructor import construct_ooh_radical
from structure_exporter import export_structure_to_cif
from pymatgen.io.ase import AseAtomsAdaptor


def test_export_structure_creates_file():
    """Test that export creates a CIF file"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    output_path = Path(__file__).parent / "test_output.cif"

    # Clean up if exists
    if output_path.exists():
        output_path.unlink()

    export_structure_to_cif(modified_atoms, output_path)

    assert output_path.exists()

    # Clean up
    output_path.unlink()


def test_export_structure_readable():
    """Test that exported CIF can be loaded back"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    output_path = Path(__file__).parent / "test_output.cif"
    export_structure_to_cif(modified_atoms, output_path)

    # Load back and verify
    loaded_structure = load_structure(output_path)
    assert len(loaded_structure) == original_num_atoms + 2

    # Clean up
    output_path.unlink()


def test_export_structure_contains_h():
    """Test that exported structure contains H atom"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    output_path = Path(__file__).parent / "test_output.cif"
    export_structure_to_cif(modified_atoms, output_path)

    loaded_structure = load_structure(output_path)
    elements = [str(site.specie) for site in loaded_structure]
    assert 'H' in elements

    # Clean up
    output_path.unlink()
