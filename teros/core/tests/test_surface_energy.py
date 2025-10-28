"""Tests for surface energy calculator module."""

import pytest
import numpy as np
from aiida import load_profile
from aiida.orm import Float, Int, StructureData
from ase import Atoms


@pytest.fixture(scope="module", autouse=True)
def load_aiida_profile():
    """Load AiiDA profile for all tests."""
    load_profile('psteros')


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
    rng = np.random.default_rng(seed=42)
    atoms = Atoms(
        symbols='Ag12P4O16',
        positions=rng.random((32, 3)) * [12.0, 8.0, 20.0],
        cell=[12.0, 8.0, 20.0],
        pbc=True
    )
    return StructureData(ase=atoms)


@pytest.fixture
def ag3po4_slab_2oh():
    """Ag3PO4 slab + 2 OH groups (Ag12P4O18H4).

    With H_{2x} formula: x=2 OH groups requires 2x=4 H atoms.
    Formula: Ag12P4O16 + 2OH -> Ag12P4O18H4
    """
    rng = np.random.default_rng(seed=42)
    atoms = Atoms(
        symbols='Ag12P4O18H4',
        positions=rng.random((38, 3)) * [12.0, 8.0, 20.0],
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


def test_get_formula_dict(ag3po4_bulk):
    """Test formula dictionary extraction from structure."""
    from teros.core.surface_hydroxylation.surface_energy import _get_formula_dict

    formula = _get_formula_dict(ag3po4_bulk)

    assert formula == {'Ag': 3, 'P': 1, 'O': 4}


def test_analyze_composition_pristine(ag3po4_bulk, ag3po4_slab_pristine):
    """Test composition analysis on pristine slab (no OH, no vacancies)."""
    from teros.core.surface_hydroxylation.surface_energy import analyze_composition

    comp = analyze_composition(ag3po4_slab_pristine, ag3po4_bulk)

    assert comp['n'] == 4
    assert comp['x'] == 0
    assert comp['y'] == 0
    assert comp['n_h'] == 0
    assert comp['n_o_deficit'] == 0


def test_analyze_composition_2oh(ag3po4_bulk, ag3po4_slab_2oh):
    """Test composition analysis on 2 OH structure.

    With H_{2x} formula: 4 H atoms → x=2 OH groups (4 / 2 = 2).
    Formula: Ag12P4O18H4 means n=4, x=2, y=0.
    """
    from teros.core.surface_hydroxylation.surface_energy import analyze_composition

    comp = analyze_composition(ag3po4_slab_2oh, ag3po4_bulk)

    assert comp['n'] == 4
    assert comp['x'] == 2  # 4 H atoms / 2 = 2 OH groups (H_{2x} formula)
    assert comp['y'] == 0  # y = (-2 + 2) / 2 = 0
    assert comp['n_h'] == 4
    assert comp['n_o_deficit'] == -2  # Net gain of 2 O from OH


def test_calc_delta_g_reaction1(expected_values):
    """Test Reaction 1 formation energy: H2O/O2 reservoirs."""
    from teros.core.surface_hydroxylation.surface_energy import calc_delta_g_reaction1

    delta_g = calc_delta_g_reaction1(
        E_slab=Float(expected_values['E_modified']),
        E_bulk=Float(expected_values['E_bulk']),
        n=Int(expected_values['n']),
        x=Int(expected_values['x']),
        y=Int(expected_values['y']),
        delta_mu_h2o=Float(expected_values['delta_mu_h2o']),
        delta_mu_o2=Float(expected_values['delta_mu_o2']),
    )

    # Should match worked example from Section S3
    # ΔG = E_slab + 2y·μ_O - n·E_bulk - x·μ_H2O
    # For pristine: E_pristine - n*E_bulk = -1600.8 - 4*(-400.2) = 0
    # For modified: -1498.5 + 0 - 4*(-400.2) - 2*(-0.5774) = 103.4548 eV
    assert abs(delta_g.value - expected_values['delta_g']) < 0.001
