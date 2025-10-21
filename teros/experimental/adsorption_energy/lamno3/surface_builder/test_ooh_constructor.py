# test_ooh_constructor.py
import pytest
import numpy as np
from pathlib import Path
from structure_loader import load_structure
from surface_detector import find_surface_mn
from coordination_finder import find_coordinated_oxygens
from oxygen_selector import select_most_exposed_oxygen
from ooh_constructor import construct_ooh_radical


def test_construct_ooh_adds_two_atoms():
    """Test that OOH construction adds exactly 2 atoms (O and H)"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    original_num_atoms = len(structure)

    # Convert to ASE for modification
    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    assert len(modified_atoms) == original_num_atoms + 2


def test_construct_ooh_adds_correct_species():
    """Test that added atoms are O and H"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    # Check that new atoms are O and H
    new_atom_symbols = modified_atoms.get_chemical_symbols()[original_num_atoms:]
    assert 'O' in new_atom_symbols
    assert 'H' in new_atom_symbols


def test_construct_ooh_bond_lengths():
    """Test that O-O and O-H bond lengths are reasonable"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    positions = modified_atoms.get_positions()

    # Get positions of base O, new O, and H
    base_o_pos = exposed_o_pos
    new_o_pos = positions[original_num_atoms]  # First new atom (O)
    h_pos = positions[original_num_atoms + 1]  # Second new atom (H)

    # Check O-O bond length (~1.45 Å for OOH)
    oo_distance = np.linalg.norm(new_o_pos - base_o_pos)
    assert 1.3 < oo_distance < 1.6, f"O-O bond length {oo_distance:.2f} Å out of range"

    # Check O-H bond length (~0.97 Å for O-H)
    oh_distance = np.linalg.norm(h_pos - new_o_pos)
    assert 0.85 < oh_distance < 1.1, f"O-H bond length {oh_distance:.2f} Å out of range"


def test_construct_ooh_orientation():
    """Test that OOH points away from surface (positive z direction)"""
    cif_path = Path(__file__).parent / "LaMnO3_100_A4_surface.cif"
    structure = load_structure(cif_path)
    mn_index, _ = find_surface_mn(structure)
    o_indices = find_coordinated_oxygens(structure, mn_index)
    exposed_o_idx, exposed_o_pos = select_most_exposed_oxygen(structure, o_indices)

    from pymatgen.io.ase import AseAtomsAdaptor
    ase_atoms = AseAtomsAdaptor.get_atoms(structure)
    original_num_atoms = len(ase_atoms)

    modified_atoms = construct_ooh_radical(ase_atoms, exposed_o_idx, exposed_o_pos)

    positions = modified_atoms.get_positions()
    base_o_z = exposed_o_pos[2]
    new_o_z = positions[original_num_atoms][2]
    h_z = positions[original_num_atoms + 1][2]

    # New O and H should be above base O
    assert new_o_z > base_o_z
    assert h_z > new_o_z
