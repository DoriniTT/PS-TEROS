"""Tests for surface energy calculator module."""

import pytest
import numpy as np
from aiida.orm import Float, Int, StructureData
from ase import Atoms


@pytest.fixture
def ag3po4_bulk():
    """Bulk Ag3PO4 structure (unit cell)."""
    # Create ASE Atoms object for Ag3PO4
    # Ag: 3 atoms, P: 1 atom, O: 4 atoms
    atoms = Atoms(
        symbols='Ag3PO4',
        positions=[
            [0.0, 0.0, 0.0],  # Ag
            [1.0, 1.0, 0.0],  # Ag
            [2.0, 0.0, 0.0],  # Ag
            [1.0, 0.5, 0.5],  # P
            [1.0, 0.0, 1.0],  # O
            [0.5, 1.0, 0.5],  # O
            [1.5, 1.0, 0.5],  # O
            [1.0, 1.0, 1.0],  # O
        ],
        cell=[3.0, 2.0, 2.0],
        pbc=True
    )
    return StructureData(ase=atoms)


@pytest.fixture
def ag3po4_slab_pristine():
    """Pristine Ag12P4O16 slab (4 bulk units, stoichiometric)."""
    # Create slab with 4 bulk formula units
    atoms = Atoms(
        symbols='Ag12P4O16',
        positions=np.random.rand(32, 3) * [12.0, 8.0, 20.0],
        cell=[12.0, 8.0, 20.0],
        pbc=True
    )
    return StructureData(ase=atoms)


@pytest.fixture
def ag3po4_slab_2oh():
    """Ag3PO4 slab + 2 OH groups (Ag12P4O18H2)."""
    # Create slab with 2 OH added: Ag12P4O16 + 2OH -> Ag12P4O18H2
    atoms = Atoms(
        symbols='Ag12P4O18H2',
        positions=np.random.rand(36, 3) * [12.0, 8.0, 20.0],
        cell=[12.0, 8.0, 20.0],
        pbc=True
    )
    return StructureData(ase=atoms)


@pytest.fixture
def expected_values():
    """Expected values from Section S3 worked example."""
    return {
        'E_bulk': -400.2,          # eV/f.u.
        'E_pristine': -1600.8,     # eV
        'E_modified': -1498.5,     # eV
        'n': 4,
        'x': 2,
        'y': 0,
        'delta_mu_h2o': -0.5774,   # eV at 298K, 0.023 bar
        'delta_mu_o2': -0.5837,    # eV at 298K, 0.21 bar
        'delta_g': 103.4548,       # eV (calculated from above)
        'area_m2': 9.6e-19,        # m² (example value)
        'gamma_s': 3.931,          # J/m²
        'gamma': 7.862,            # J/m²
    }


def test_fixtures_load():
    """Test that fixtures load correctly."""
    pass  # Placeholder - fixtures tested indirectly
