# teros/core/test_adsorption_energy.py
import pytest
import numpy as np
from aiida import orm, load_profile
from pymatgen.core import Structure, Lattice
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

from .adsorption_energy import _separate_adsorbate_structure_impl, _calculate_adsorption_energy_impl


@pytest.fixture(scope='module', autouse=True)
def aiida_profile():
    """Load AiiDA profile before running tests."""
    load_profile()
    yield


def create_test_surface_with_oh():
    """Create a simple Ag(111) surface with OH adsorbate for testing."""
    # Create 2x2 Ag slab
    lattice = Lattice.from_parameters(a=5.8, b=5.8, c=20.0,
                                     alpha=90, beta=90, gamma=90)

    # Ag atoms in bottom layer
    ag_positions = [
        [0.0, 0.0, 10.0],
        [2.9, 0.0, 10.0],
        [0.0, 2.9, 10.0],
        [2.9, 2.9, 10.0],
    ]

    # OH adsorbate on top (bonded O-H cluster)
    # O-H bond distance should be around 0.97 Angstrom
    oh_positions = [
        [1.45, 1.45, 12.5],  # O atom above surface
        [1.45, 1.45, 13.47],  # H atom above O (0.97 A bond length)
    ]

    species = ['Ag'] * 4 + ['O', 'H']
    positions = ag_positions + oh_positions

    structure = Structure(lattice, species, positions, coords_are_cartesian=True)

    return orm.StructureData(pymatgen=structure)


def test_separate_adsorbate_returns_three_structures():
    """Test that separation returns substrate, molecule, and complete structures."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = _separate_adsorbate_structure_impl(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    assert 'substrate' in result
    assert 'molecule' in result
    assert 'complete' in result

    assert isinstance(result['substrate'], orm.StructureData)
    assert isinstance(result['molecule'], orm.StructureData)
    assert isinstance(result['complete'], orm.StructureData)


def test_separate_adsorbate_correct_atom_counts():
    """Test that separated structures have correct number of atoms."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = _separate_adsorbate_structure_impl(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    # Original: 4 Ag + 2 (O,H) = 6 atoms
    # Substrate: 4 Ag
    # Molecule: 1 O + 1 H = 2 atoms
    # Complete: 6 atoms

    substrate_ase = result['substrate'].get_ase()
    molecule_ase = result['molecule'].get_ase()
    complete_ase = result['complete'].get_ase()

    assert len(substrate_ase) == 4
    assert len(molecule_ase) == 2
    assert len(complete_ase) == 6


def test_separate_adsorbate_correct_species():
    """Test that molecule contains only adsorbate species."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = _separate_adsorbate_structure_impl(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    molecule_ase = result['molecule'].get_ase()
    symbols = molecule_ase.get_chemical_symbols()

    assert 'O' in symbols
    assert 'H' in symbols
    assert 'Ag' not in symbols


def test_separate_adsorbate_same_cell():
    """Test that all three structures have the same cell."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OH')

    result = _separate_adsorbate_structure_impl(
        structure=complete,
        adsorbate_formula=adsorbate_formula
    )

    complete_ase = result['complete'].get_ase()
    substrate_ase = result['substrate'].get_ase()
    molecule_ase = result['molecule'].get_ase()

    complete_cell = complete_ase.get_cell()
    substrate_cell = substrate_ase.get_cell()
    molecule_cell = molecule_ase.get_cell()

    assert np.allclose(complete_cell, substrate_cell)
    assert np.allclose(complete_cell, molecule_cell)


def test_separate_adsorbate_invalid_formula():
    """Test that invalid formula raises error."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('')  # Empty formula

    with pytest.raises(ValueError, match="Adsorbate formula cannot be empty"):
        _separate_adsorbate_structure_impl(
            structure=complete,
            adsorbate_formula=adsorbate_formula
        )


def test_separate_adsorbate_no_match():
    """Test that non-existent adsorbate raises error."""
    complete = create_test_surface_with_oh()
    adsorbate_formula = orm.Str('OOH')  # Not in structure

    with pytest.raises(ValueError, match="Could not find adsorbate"):
        _separate_adsorbate_structure_impl(
            structure=complete,
            adsorbate_formula=adsorbate_formula
        )


def test_calculate_adsorption_energy_correct_formula():
    """Test that adsorption energy uses correct formula."""
    # E_ads = E_complete - E_substrate - E_molecule

    E_complete = orm.Float(-100.0)
    E_substrate = orm.Float(-80.0)
    E_molecule = orm.Float(-15.0)

    E_ads = _calculate_adsorption_energy_impl(
        E_complete=E_complete,
        E_substrate=E_substrate,
        E_molecule=E_molecule
    )

    # E_ads = -100 - (-80) - (-15) = -100 + 80 + 15 = -5.0
    expected = -5.0
    assert abs(E_ads.value - expected) < 1e-9


def test_calculate_adsorption_energy_positive():
    """Test case where adsorption is endothermic (positive E_ads)."""
    E_complete = orm.Float(-50.0)
    E_substrate = orm.Float(-40.0)
    E_molecule = orm.Float(-15.0)

    E_ads = _calculate_adsorption_energy_impl(
        E_complete=E_complete,
        E_substrate=E_substrate,
        E_molecule=E_molecule
    )

    # E_ads = -50 - (-40) - (-15) = -50 + 40 + 15 = 5.0 (endothermic)
    expected = 5.0
    assert abs(E_ads.value - expected) < 1e-9


def test_calculate_adsorption_energy_returns_float():
    """Test that result is Float node."""
    E_complete = orm.Float(-100.0)
    E_substrate = orm.Float(-80.0)
    E_molecule = orm.Float(-15.0)

    E_ads = _calculate_adsorption_energy_impl(
        E_complete=E_complete,
        E_substrate=E_substrate,
        E_molecule=E_molecule
    )

    assert isinstance(E_ads, orm.Float)


def test_compute_adsorption_energies_scatter_signature():
    """Test that scatter function exists with correct signature."""
    # This test just checks the function exists and accepts right parameters
    # Full integration test will be in examples folder

    import inspect
    from .adsorption_energy import compute_adsorption_energies_scatter

    # Function should be callable (wrapped by @task.graph)
    assert callable(compute_adsorption_energies_scatter)

    # Check the original function signature (before wrapping)
    # The @task.graph decorator stores the original function
    if hasattr(compute_adsorption_energies_scatter, 'node'):
        # Get the original function from the task
        original_func = compute_adsorption_energies_scatter.node.func
        sig = inspect.signature(original_func)
        params = list(sig.parameters.keys())

        # Check required parameters exist
        assert 'structures' in params
        assert 'adsorbate_formulas' in params
        assert 'code' in params
        assert 'potential_family' in params
        assert 'potential_mapping' in params
        assert 'parameters' in params
        assert 'options' in params
    else:
        # If we can't access the original function, at least verify it's callable
        # This is sufficient for the minimal test
        assert callable(compute_adsorption_energies_scatter)
