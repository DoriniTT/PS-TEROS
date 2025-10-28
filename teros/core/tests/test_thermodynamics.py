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
