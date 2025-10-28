"""Test that thermodynamics module can be imported easily."""


def test_import_from_surface_hydroxylation():
    """Test import from main module."""
    from teros.core.surface_hydroxylation import JanafDatabase

    db = JanafDatabase()
    assert db is not None


def test_direct_import():
    """Test direct import from thermodynamics."""
    from teros.core.surface_hydroxylation.thermodynamics import (
        JanafDatabase,
        KJ_TO_EV,
        KB_EV,
    )

    assert JanafDatabase is not None
    assert KJ_TO_EV > 0
    assert KB_EV > 0
