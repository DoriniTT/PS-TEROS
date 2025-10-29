"""Tests for WorkGraph integration of surface energy calculations."""
import pytest
import numpy as np
from unittest.mock import Mock
from aiida import orm, load_profile


@pytest.fixture(scope="module", autouse=True)
def load_aiida_profile():
    """Load AiiDA profile for all tests."""
    load_profile('psteros')


def test_create_surface_energy_task():
    """Test that surface energy task is created with correct structure."""
    from teros.core.surface_hydroxylation.surface_energy_workgraph import (
        create_surface_energy_task
    )

    # Mock upstream tasks
    aggregate_task = Mock()
    aggregate_task.outputs = Mock()
    aggregate_task.outputs.structures = Mock()
    aggregate_task.outputs.energies = Mock()
    bulk_task = Mock()
    bulk_task.outputs = Mock()
    bulk_task.outputs.structure = Mock()
    bulk_task.outputs.energy = Mock()
    pristine_task = Mock()
    pristine_task.outputs = Mock()

    # Create task (returns TaskSocketNamespace representing outputs)
    task_outputs = create_surface_energy_task(
        aggregate_results_task=aggregate_task,
        bulk_relax_task=bulk_task,
        pristine_relax_task=pristine_task,
    )

    # Verify task outputs structure (WorkGraph tasks return outputs namespace directly)
    assert task_outputs is not None
    assert hasattr(task_outputs, 'reaction1_results')
    assert hasattr(task_outputs, 'reaction2_results')
    assert hasattr(task_outputs, 'reaction3_results')


def test_calculate_all_surface_energies_with_fixtures():
    """Test that task executes correctly with real fixtures."""
    from teros.core.surface_hydroxylation.surface_energy_workgraph import (
        _calculate_all_surface_energies_impl
    )
    from aiida import orm

    # Create test structures with proper Ag3PO4 stoichiometry
    # Use the same fixtures as test_surface_energy.py for consistency

    # Pristine: Ag12P4O16 slab (4 bulk formula units, stoichiometric)
    from ase import Atoms
    rng = np.random.default_rng(seed=42)
    pristine_atoms = Atoms(
        symbols='Ag12P4O16',
        positions=rng.random((32, 3)) * [12.0, 8.0, 20.0],
        cell=[12.0, 8.0, 20.0],
        pbc=True
    )
    pristine_structure = orm.StructureData(ase=pristine_atoms)

    # Modified: Ag12P4O18H4 (2 OH groups added, x=2)
    # Formula: Ag12P4O16 + 2OH -> Ag12P4O18H4
    modified_atoms = Atoms(
        symbols='Ag12P4O18H4',
        positions=rng.random((38, 3)) * [12.0, 8.0, 20.0],
        cell=[12.0, 8.0, 20.0],
        pbc=True
    )
    modified_structure = orm.StructureData(ase=modified_atoms)

    # Bulk structure: Ag3PO4 unit cell
    bulk_atoms = Atoms(
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
    bulk_structure = orm.StructureData(ase=bulk_atoms)

    # Test data
    structures_dict = {
        'pristine': pristine_structure,
        'modified_2oh': modified_structure,
    }
    energies_dict = {
        'pristine': orm.Float(-1600.8),
        'modified_2oh': orm.Float(-1498.5),
    }
    bulk_energy = orm.Float(-400.2)
    temperature = 298.0
    pressures = {'H2O': 0.023, 'O2': 0.21, 'H2': 1.0}

    # Execute task implementation function
    results = _calculate_all_surface_energies_impl(
        structures_dict=structures_dict,
        energies_dict=energies_dict,
        bulk_structure=bulk_structure,
        bulk_energy=bulk_energy,
        temperature=temperature,
        pressures=pressures,
    )

    # Verify all reactions completed
    assert 'reaction1_results' in results
    assert 'reaction2_results' in results
    assert 'reaction3_results' in results

    # When calling the task function directly (not in WorkGraph),
    # it returns the dict values directly from the function body
    r1 = results['reaction1_results']
    r2 = results['reaction2_results']
    r3 = results['reaction3_results']

    # If they are AiidaDict objects, extract the dict
    if hasattr(r1, 'get_dict'):
        r1 = r1.get_dict()
        r2 = r2.get_dict()
        r3 = r3.get_dict()

    # Verify structure of each reaction result
    for reaction, result in [('reaction1', r1), ('reaction2', r2), ('reaction3', r3)]:
        assert 'surface_energies' in result, f"{reaction} missing surface_energies"
        assert 'formation_energies' in result, f"{reaction} missing formation_energies"
        assert 'stoichiometry' in result, f"{reaction} missing stoichiometry"
        assert 'reference_data' in result, f"{reaction} missing reference_data"
        assert 'error' not in result, f"{reaction} has error: {result.get('error')}"

    # Verify pristine is included
    assert 'pristine' in r1['surface_energies']
    assert 'modified_2oh' in r1['surface_energies']

    # Verify surface energies are non-negative (J/m²)
    # Pristine is reference, so γ₀ = 0
    assert r1['surface_energies']['pristine'] >= 0
    # Modified structures should have calculated surface energies
    assert isinstance(r1['surface_energies']['modified_2oh'], (int, float))
