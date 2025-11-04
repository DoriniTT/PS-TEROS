"""Tests for AIMD functions."""
import pytest
from teros.core.aimd_functions import prepare_aimd_parameters


def test_prepare_aimd_parameters_required_only():
    """Test prepare_aimd_parameters with only required TEBEG and NSW."""
    base_params = {
        'PREC': 'Normal',
        'ENCUT': 400,
        'EDIFF': 1e-5,
    }

    stage_config = {
        'TEBEG': 300,
        'NSW': 100,
    }

    result = prepare_aimd_parameters(base_params, stage_config)

    # Should have base parameters
    assert result['PREC'] == 'Normal'
    assert result['ENCUT'] == 400
    assert result['EDIFF'] == 1e-5

    # Should have AIMD parameters from stage
    assert result['TEBEG'] == 300
    assert result['TEEND'] == 300  # Defaults to TEBEG
    assert result['NSW'] == 100
    assert result['IBRION'] == 0  # Forced to 0 for MD


def test_prepare_aimd_parameters_with_optional():
    """Test prepare_aimd_parameters with optional AIMD parameters."""
    base_params = {
        'PREC': 'Normal',
        'POTIM': 1.0,  # Base timestep
    }

    stage_config = {
        'TEBEG': 300,
        'NSW': 100,
        'TEEND': 400,  # Different end temp
        'POTIM': 2.0,  # Override base
        'MDALGO': 2,
        'SMASS': 0.0,
    }

    result = prepare_aimd_parameters(base_params, stage_config)

    # Stage AIMD params should override base
    assert result['TEBEG'] == 300
    assert result['TEEND'] == 400
    assert result['NSW'] == 100
    assert result['POTIM'] == 2.0  # Overridden
    assert result['MDALGO'] == 2
    assert result['SMASS'] == 0.0


def test_prepare_aimd_parameters_missing_tebeg():
    """Test that ValueError is raised if TEBEG missing."""
    base_params = {'PREC': 'Normal'}
    stage_config = {'NSW': 100}

    with pytest.raises(ValueError, match="must contain 'TEBEG' and 'NSW'"):
        prepare_aimd_parameters(base_params, stage_config)


def test_prepare_aimd_parameters_missing_nsw():
    """Test that ValueError is raised if NSW missing."""
    base_params = {'PREC': 'Normal'}
    stage_config = {'TEBEG': 300}

    with pytest.raises(ValueError, match="must contain 'TEBEG' and 'NSW'"):
        prepare_aimd_parameters(base_params, stage_config)
