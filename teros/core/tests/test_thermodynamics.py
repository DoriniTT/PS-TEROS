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
