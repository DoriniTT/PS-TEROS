"""Tests for JANAF thermodynamics database."""

import pytest
from pathlib import Path
from teros.core.surface_hydroxylation.thermodynamics import JanafDatabase


def test_database_loads():
    """Test that database loads without errors."""
    db = JanafDatabase()
    assert db is not None


def test_available_species():
    """Test available species list."""
    db = JanafDatabase()
    species = db.available_species

    assert 'H2O' in species
    assert 'H2' in species
    assert 'O2' in species
    assert len(species) == 3


def test_get_mu_correction_basic():
    """Test basic Δμ⁰ lookup at 298 K."""
    db = JanafDatabase()

    # Get correction at 298 K
    mu_h2o = db.get_mu_correction('H2O', T=298)
    mu_h2 = db.get_mu_correction('H2', T=298)
    mu_o2 = db.get_mu_correction('O2', T=298)

    # All should be negative (favorable)
    assert mu_h2o < 0
    assert mu_h2 < 0
    assert mu_o2 < 0

    # Should be within reasonable range (-1 to 0 eV)
    assert -1.0 < mu_h2o < 0
    assert -1.0 < mu_h2 < 0
    assert -1.0 < mu_o2 < 0


def test_get_mu_correction_temperature_dependence():
    """Test that Δμ⁰ becomes more negative with temperature."""
    db = JanafDatabase()

    mu_298 = db.get_mu_correction('H2O', T=298)
    mu_500 = db.get_mu_correction('H2O', T=500)

    # Higher temperature should give more negative correction
    assert mu_500 < mu_298


def test_invalid_species():
    """Test error on invalid species."""
    db = JanafDatabase()

    with pytest.raises(ValueError, match="not found"):
        db.get_mu_correction('CO2', T=298)


def test_invalid_temperature():
    """Test error on invalid temperature."""
    db = JanafDatabase()

    # Not a 50 K step
    with pytest.raises(ValueError, match="not found in database"):
        db.get_mu_correction('H2O', T=299)

    # Out of range
    with pytest.raises(ValueError, match="not found in database"):
        db.get_mu_correction('H2O', T=3000)


def test_temperature_must_be_discrete():
    """Test that intermediate temperatures are not allowed."""
    db = JanafDatabase()

    # 325 K is between 300 and 350
    with pytest.raises(ValueError, match="not found in database"):
        db.get_mu_correction('H2O', T=325)


def test_get_raw_data():
    """Test raw JANAF data retrieval."""
    db = JanafDatabase()

    raw = db.get_raw_data('H2O', T=298)

    assert 'H' in raw
    assert 'S' in raw
    assert 'delta_mu' in raw

    # Enthalpy should be positive and reasonable
    assert 0 < raw['H'] < 20  # kJ/mol

    # Entropy should be positive
    assert raw['S'] > 0  # J/(mol·K)

    # Delta mu should match correction method
    mu = db.get_mu_correction('H2O', T=298)
    assert abs(raw['delta_mu'] - mu) < 1e-6


def test_list_temperatures():
    """Test temperature listing."""
    db = JanafDatabase()

    temps = db.list_temperatures()

    # Should have multiple points covering 0-2000 K range
    assert len(temps) > 20  # At least 20+ points

    # Should be sorted
    assert temps == sorted(temps)

    # Check range
    assert temps[0] == 0
    assert temps[-1] == 2000

    # Check specific important temperatures exist
    assert 50 in temps
    assert 298 in temps  # Room temperature
    assert 300 in temps
    assert 500 in temps
    assert 1000 in temps
