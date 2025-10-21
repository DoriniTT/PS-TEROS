"""Tests for surface_hydroxylation utilities."""

import pytest
from aiida import load_profile
from aiida.orm import StructureData
from ase import Atoms
from teros.core.surface_hydroxylation.utils import aiida_to_ase, ase_to_aiida

# Load AiiDA profile before tests
load_profile('psteros')


def test_aiida_to_ase_conversion():
    """Test StructureData → ASE Atoms conversion."""
    # Create simple AiiDA structure (2-atom H2 molecule)
    structure = StructureData(cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    structure.append_atom(position=(5.0, 5.0, 5.0), symbols='H')
    structure.append_atom(position=(5.5, 5.0, 5.0), symbols='H')

    # Convert to ASE
    atoms = aiida_to_ase(structure)

    # Verify
    assert isinstance(atoms, Atoms)
    assert len(atoms) == 2
    assert atoms.get_chemical_symbols() == ['H', 'H']


def test_ase_to_aiida_conversion():
    """Test ASE Atoms → StructureData conversion."""
    # Create ASE Atoms
    atoms = Atoms('H2', positions=[[0, 0, 0], [0.74, 0, 0]], cell=[10, 10, 10], pbc=True)

    # Convert to AiiDA
    structure = ase_to_aiida(atoms)

    # Verify
    assert isinstance(structure, StructureData)
    assert len(structure.sites) == 2
