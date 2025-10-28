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
        'gamma_s': 8.633,          # J/m² (gamma_s = ΔG/(2A) = 103.4548*1.602e-19/(2*9.6e-19))
        'gamma': 17.266,           # J/m² (gamma = 2*gamma_s - gamma_0, where gamma_0=0 for pristine)
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


def test_calc_delta_g_reaction2():
    """Test Reaction 2 formation energy: H2/H2O reservoirs."""
    from teros.core.surface_hydroxylation.surface_energy import calc_delta_g_reaction2

    # Use simple values for validation
    delta_g = calc_delta_g_reaction2(
        E_slab=Float(-1500.0),
        E_bulk=Float(-400.0),
        n=Int(4),
        x=Int(2),
        y=Int(0),
        delta_mu_h2=Float(-0.3),
        delta_mu_h2o=Float(-0.5),
    )

    # ΔG = E + (2y-x)·μ_H2O - n·E_bulk - x·μ_H2
    # = -1500 + (0-2)*(-0.5) - 4*(-400) - 2*(-0.3)
    # = -1500 + 1.0 + 1600 + 0.6 = 101.6
    expected = -1500.0 + (-2)*(-0.5) - 4*(-400.0) - 2*(-0.3)
    assert abs(delta_g.value - expected) < 0.001


def test_calc_delta_g_reaction3():
    """Test Reaction 3 formation energy: H2/O2 reservoirs."""
    from teros.core.surface_hydroxylation.surface_energy import calc_delta_g_reaction3

    # Use simple values for validation
    delta_g = calc_delta_g_reaction3(
        E_slab=Float(-1500.0),
        E_bulk=Float(-400.0),
        n=Int(4),
        x=Int(2),
        y=Int(0),
        delta_mu_h2=Float(-0.3),
        delta_mu_o2=Float(-0.6),
    )

    # ΔG = E + (2y-x)·μ_O - n·E_bulk - x·μ_H2
    # μ_O = 0.5·Δμ(O2) = 0.5*(-0.6) = -0.3
    # = -1500 + (-2)*(-0.3) - 4*(-400) - 2*(-0.3)
    # = -1500 + 0.6 + 1600 + 0.6 = 101.2
    mu_o = 0.5 * (-0.6)
    expected = -1500.0 + (-2)*mu_o - 4*(-400.0) - 2*(-0.3)
    assert abs(delta_g.value - expected) < 0.001


def test_select_reaction_function():
    """Test reaction function selector."""
    from teros.core.surface_hydroxylation.surface_energy import (
        select_reaction_function,
        calc_delta_g_reaction1,
        calc_delta_g_reaction2,
        calc_delta_g_reaction3,
    )

    assert select_reaction_function(1) == calc_delta_g_reaction1
    assert select_reaction_function(2) == calc_delta_g_reaction2
    assert select_reaction_function(3) == calc_delta_g_reaction3

    # Test invalid reaction number
    with pytest.raises(ValueError, match="must be 1, 2, or 3"):
        select_reaction_function(4)


def test_get_surface_area(ag3po4_slab_pristine):
    """Test surface area calculation from cell parameters."""
    from teros.core.surface_hydroxylation.surface_energy import get_surface_area

    area = get_surface_area(ag3po4_slab_pristine)

    # Area should be cross product magnitude of a and b vectors
    # Cell is [12.0, 8.0, 20.0], so area = 12*8 = 96 Ų = 96e-20 m²
    assert isinstance(area, Float)
    assert area.value > 0
    assert abs(area.value - 96e-20) < 1e-21  # Tolerance for float comparison


def test_calc_gamma_s(expected_values):
    """Test average surface energy calculation (Equation 9)."""
    from teros.core.surface_hydroxylation.surface_energy import calc_gamma_s

    # γ_s = ΔG / (2A)
    # Using expected values from Section S3
    gamma_s = calc_gamma_s(
        delta_g=Float(expected_values['delta_g']),
        area=Float(expected_values['area_m2'])
    )

    assert isinstance(gamma_s, Float)
    assert abs(gamma_s.value - expected_values['gamma_s']) < 0.01


def test_calc_gamma(expected_values):
    """Test corrected surface energy calculation (Equation 10)."""
    from teros.core.surface_hydroxylation.surface_energy import calc_gamma

    # γ = 2γ_s - γ_0
    # For single-sided modification, removes pristine surface contribution
    gamma = calc_gamma(
        gamma_s_modified=Float(expected_values['gamma_s']),
        gamma_0_pristine=Float(0.0)  # Assuming pristine is reference
    )

    assert isinstance(gamma, Float)
    assert abs(gamma.value - expected_values['gamma']) < 0.01


def test_calc_gamma_with_pristine():
    """Test gamma correction with non-zero pristine."""
    from teros.core.surface_hydroxylation.surface_energy import calc_gamma

    # If pristine has γ_0 = 1.0 J/m² and modified has γ_s = 3.0 J/m²
    # Then γ = 2*3.0 - 1.0 = 5.0 J/m²
    gamma = calc_gamma(
        gamma_s_modified=Float(3.0),
        gamma_0_pristine=Float(1.0)
    )

    assert abs(gamma.value - 5.0) < 1e-6


def test_calculate_surface_energies_full_workflow(
    ag3po4_bulk,
    ag3po4_slab_pristine,
    ag3po4_slab_2oh,
    expected_values
):
    """Test complete orchestration function with pristine + modified structures."""
    from teros.core.surface_hydroxylation.surface_energy import calculate_surface_energies

    structures = {
        'pristine': ag3po4_slab_pristine,
        '2oh': ag3po4_slab_2oh,
    }
    energies = {
        'pristine': Float(expected_values['E_pristine']),
        '2oh': Float(expected_values['E_modified']),
    }

    results = calculate_surface_energies(
        structures_dict=structures,
        energies_dict=energies,
        bulk_structure=ag3po4_bulk,
        bulk_energy=Float(expected_values['E_bulk']),
        temperature=298.0,
        pressures={'H2O': 0.023, 'O2': 0.21},
        which_reaction=1,
    )

    # Check structure of results
    assert 'surface_energies' in results
    assert 'formation_energies' in results
    assert 'stoichiometry' in results
    assert 'reference_data' in results

    # Check stoichiometry was calculated correctly
    assert results['stoichiometry']['2oh']['n'] == 4
    assert results['stoichiometry']['2oh']['x'] == 2
    assert results['stoichiometry']['2oh']['y'] == 0

    # Check that pristine was identified
    assert results['stoichiometry']['pristine']['x'] == 0
    assert results['stoichiometry']['pristine']['y'] == 0

    # Check reference data
    assert results['reference_data']['temperature'] == 298.0
    assert results['reference_data']['reaction_used'] == 1
    assert 'gamma_0' in results['reference_data']
