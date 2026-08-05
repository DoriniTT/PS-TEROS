"""
Pytest configuration and fixtures for PS-TEROS tests.
"""

import pytest
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# =============================================================================
# MARKERS
# =============================================================================

def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "requires_aiida: marks tests that require AiiDA to be configured"
    )
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "tier1: fast unit tests without AiiDA"
    )
    config.addinivalue_line(
        "markers", "tier2: integration tests with AiiDA mock database"
    )
    config.addinivalue_line(
        "markers", "tier3: end-to-end tests with real calculations"
    )
    config.addinivalue_line(
        "markers", "localwork: tier3 tests using VASP-6.5.1@localwork"
    )
    config.addinivalue_line(
        "markers", "obelix: tier3 tests using VASP-6.5.1-idefix-4@obelix"
    )


# =============================================================================
# FIXTURES
# =============================================================================

@pytest.fixture
def sample_formula_strings():
    """Sample chemical formulas for testing."""
    return {
        'simple': ['H2O', 'O2', 'H2', 'CO2', 'NH3'],
        'complex': ['OOH', 'OH', 'H2O2', 'CH3OH'],
        'oxides': ['Ag2O', 'Fe2O3', 'TiO2', 'Al2O3'],
        'ternary': ['Ag3PO4', 'LaMnO3', 'BaTiO3'],
    }


@pytest.fixture
def sample_dict_for_merge():
    """Sample dictionaries for merge testing."""
    return {
        'base': {
            'a': 1,
            'b': {'c': 2, 'd': 3},
            'e': [1, 2, 3],
        },
        'override': {
            'b': {'c': 99, 'f': 100},
            'g': 'new',
        },
        'expected': {
            'a': 1,
            'b': {'c': 99, 'd': 3, 'f': 100},
            'e': [1, 2, 3],
            'g': 'new',
        },
    }


@pytest.fixture
def workflow_presets_list():
    """List of all expected workflow presets."""
    return [
        'bulk_only',
        'formation_enthalpy_only',
        'surface_thermodynamics',
        'surface_thermodynamics_unrelaxed',
        'cleavage_only',
        'relaxation_energy_only',
        'electronic_structure_bulk_only',
        'electronic_structure_slabs_only',
        'electronic_structure_bulk_and_slabs',
        'aimd_only',
        'adsorption_energy',
        'comprehensive',
    ]


# =============================================================================
# AIIDA AVAILABILITY CHECK
# =============================================================================

def check_aiida_available():
    """Check if AiiDA is properly configured."""
    try:
        from aiida import load_profile
        load_profile()
        return True
    except Exception:
        return False


AIIDA_AVAILABLE = check_aiida_available()


@pytest.fixture
def skip_without_aiida():
    """Skip test if AiiDA is not available."""
    if not AIIDA_AVAILABLE:
        pytest.skip("AiiDA not configured")


# =============================================================================
# STRUCTURE FIXTURES (require AiiDA)
# =============================================================================

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), 'fixtures')
STRUCTURES_DIR = os.path.join(FIXTURES_DIR, 'structures')


def _load_structure_from_file(filename):
    """Load an AiiDA StructureData from a VASP file in the fixtures directory."""
    from ase.io import read as ase_read
    from aiida import orm

    filepath = os.path.join(STRUCTURES_DIR, filename)
    atoms = ase_read(filepath, format='vasp')
    return orm.StructureData(ase=atoms)


@pytest.fixture
def si_diamond_structure():
    """2-atom Si diamond StructureData (requires AiiDA)."""
    if not AIIDA_AVAILABLE:
        pytest.skip("AiiDA not configured")
    return _load_structure_from_file('si_diamond.vasp')


@pytest.fixture
def sno2_rutile_structure():
    """6-atom SnO2 rutile StructureData (requires AiiDA)."""
    if not AIIDA_AVAILABLE:
        pytest.skip("AiiDA not configured")
    return _load_structure_from_file('sno2_rutile.vasp')

