"""Tests for surface_hydroxylation relaxations module."""

import pytest
from aiida import orm, load_profile


# Load profile before tests
@pytest.fixture(scope='module', autouse=True)
def load_aiida_profile():
    """Load AiiDA profile for testing."""
    load_profile('psteros')


def test_import_relax_slabs_with_semaphore():
    """Test that relax_slabs_with_semaphore can be imported."""
    from teros.core.surface_hydroxylation.relaxations import relax_slabs_with_semaphore

    assert relax_slabs_with_semaphore is not None
    assert callable(relax_slabs_with_semaphore)


def test_import_extract_total_energy():
    """Test that extract_total_energy can be imported."""
    from teros.core.surface_hydroxylation.relaxations import extract_total_energy

    assert extract_total_energy is not None
    assert callable(extract_total_energy)


def test_import_get_exit_info():
    """Test that get_exit_info can be imported."""
    from teros.core.surface_hydroxylation.relaxations import get_exit_info

    assert get_exit_info is not None
    assert callable(get_exit_info)


def test_extract_total_energy_basic():
    """Test extract_total_energy with basic energy dict."""
    from teros.core.surface_hydroxylation.relaxations import extract_total_energy
    from aiida.engine import run

    # Create test energy dict
    energies = orm.Dict(dict={
        'total_energies': {
            'energy_extrapolated': -123.45
        }
    })

    # Run calcfunction
    result = run(extract_total_energy, energies=energies)

    # Verify result
    assert isinstance(result, orm.Float)
    assert abs(result.value - (-123.45)) < 1e-6


def test_get_exit_info_success():
    """Test get_exit_info with successful exit."""
    from teros.core.surface_hydroxylation.relaxations import get_exit_info
    from aiida.engine import run

    # Create test inputs
    exit_status = orm.Int(0)
    exit_message = orm.Str('')

    # Run calcfunction
    result = run(get_exit_info, exit_status=exit_status, exit_message=exit_message)

    # Verify result
    assert 'exit_status' in result
    assert 'error' in result
    assert result['exit_status'].value == 0
    assert result['error'].value == ''


def test_get_exit_info_failure():
    """Test get_exit_info with failed exit."""
    from teros.core.surface_hydroxylation.relaxations import get_exit_info
    from aiida.engine import run

    # Create test inputs
    exit_status = orm.Int(400)
    exit_message = orm.Str('VASP calculation did not converge')

    # Run calcfunction
    result = run(get_exit_info, exit_status=exit_status, exit_message=exit_message)

    # Verify result
    assert 'exit_status' in result
    assert 'error' in result
    assert result['exit_status'].value == 400
    assert 'converge' in result['error'].value.lower()
