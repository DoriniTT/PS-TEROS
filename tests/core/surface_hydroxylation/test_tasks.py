"""Tests for surface_hydroxylation task functions."""

import pytest
from aiida import orm, load_profile
from aiida.engine import run
from ase.build import fcc111
from teros.core.surface_hydroxylation.tasks import generate_structures


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
