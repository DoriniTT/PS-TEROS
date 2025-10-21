"""Tests for surface_hydroxylation task functions."""

import pytest
from aiida import orm, load_profile
from aiida.engine import run
from ase.build import fcc111
from teros.core.surface_hydroxylation.tasks import generate_structures, collect_results


# Load AiiDA profile once at module level
load_profile('psteros')


def test_generate_structures_hydrogen_mode():
    """Test structure generation in hydrogen mode."""
    # Create simple test slab (2x2 Pt(111) with O adlayer)
    slab = fcc111('Pt', size=(2, 2, 4), vacuum=10.0)

    # Add O atoms on top manually
    # Get the top surface atoms of Pt slab
    z_max = slab.positions[:, 2].max()
    # Add O atoms above the Pt surface
    from ase import Atoms
    o_layer = Atoms('O4',
                    positions=[[0, 0, z_max + 2.0],
                               [2.77, 0, z_max + 2.0],
                               [0, 2.77, z_max + 2.0],
                               [2.77, 2.77, z_max + 2.0]],
                    cell=slab.cell,
                    pbc=True)
    slab.extend(o_layer)
    slab.center(vacuum=10.0, axis=2)

    structure = orm.StructureData(ase=slab)

    params = orm.Dict({
        'mode': 'hydrogen',
        'species': 'O',
        'z_window': 0.5,
        'which_surface': 'top',
        'oh_dist': 0.98,
        'include_empty': False,
        'deduplicate_by_coverage': True,
        'coverage_bins': 3
    })

    # Run calcfunction
    result = run(generate_structures, structure=structure, params=params)

    # Verify outputs
    assert 'manifest' in result
    assert isinstance(result['manifest'], orm.Dict)

    # Count structure outputs (structure_0, structure_1, etc.)
    structure_keys = [k for k in result.keys() if k.startswith('structure_')]
    assert len(structure_keys) > 0, "Should have at least one structure"
    assert all(isinstance(result[k], orm.StructureData) for k in structure_keys)


def test_collect_results_with_mixed_success_failure():
    """Test result collection with some successes and failures."""
    # Create manifest
    manifest = orm.Dict({
        'variants': [
            {'name': 'oh_001_0.5', 'OH_coverage': 0.5},
            {'name': 'oh_002_1.0', 'OH_coverage': 1.0}
        ]
    })

    # Create mock structures and energies for successful relaxation
    success_structure = orm.StructureData(cell=[[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    success_structure.append_atom(position=(5.0, 5.0, 5.0), symbols='H')
    success_structure.store()
    success_energy = orm.Float(-123.45)
    success_energy.store()

    # Pass relaxation results as namespace kwargs
    # In real usage, these come from RelaxationsWorkGraph outputs
    # Format: structure_0, energy_0, exit_status_0, etc.
    relaxations_namespace = {
        'structure_0': success_structure,
        'energy_0': success_energy,
        'exit_status_0': orm.Int(0),
        'exit_status_1': orm.Int(400),
        'error_1': orm.Str('Convergence failed')
    }

    # Run collect_results with **kwargs to pass namespace
    result = run(collect_results, manifest=manifest, **relaxations_namespace)

    # Verify outputs
    assert 'successful_relaxations' in result
    assert 'failed_relaxations' in result
    assert 'statistics' in result

    # Check statistics
    stats = result['statistics'].get_dict()
    assert stats['total'] == 2
    assert stats['succeeded'] == 1
    assert stats['failed'] == 1

    # Check successful relaxations
    successful = result['successful_relaxations'].get_dict()['results']
    assert len(successful) == 1
    assert successful[0]['name'] == 'oh_001_0.5'
    assert successful[0]['energy'] == -123.45
    # Verify structure_pk is present and is an integer
    assert 'structure_pk' in successful[0]
    assert isinstance(successful[0]['structure_pk'], int)

    # Check failed relaxations
    failed = result['failed_relaxations'].get_dict()['results']
    assert len(failed) == 1
    assert failed[0]['name'] == 'oh_002_1.0'
    assert failed[0]['exit_status'] == 400
    assert failed[0]['error_message'] == 'Convergence failed'

    # Verify we can load the structure using the stored PK
    loaded_structure = orm.load_node(successful[0]['structure_pk'])
    assert isinstance(loaded_structure, orm.StructureData)
    # Verify it has the expected properties (cell and atom count)
    assert len(loaded_structure.sites) == 1
