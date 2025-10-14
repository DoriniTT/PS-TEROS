#!/usr/bin/env python
"""
Example: Selective Slab Electronic Properties Calculation

This example demonstrates how to calculate DOS and band structures
for selected slab terminations with per-slab parameter overrides.

Based on: docs/plans/2025-10-12-slab-electronic-properties.md
"""

from aiida import orm
from teros.core.workgraph import build_core_workgraph
from teros.core.builders.electronic_properties_builder import (
    get_electronic_properties_defaults,
    get_slab_electronic_properties_defaults,
)

# Get defaults from builders
bulk_defaults = get_electronic_properties_defaults()
slab_defaults = get_slab_electronic_properties_defaults()

# Custom parameters for term_2 (example: higher accuracy)
slab_custom = get_slab_electronic_properties_defaults(
    energy_cutoff=600,
    kpoints_mesh_density=0.2,  # Even denser
    band_kpoints_distance=0.1,
)

# Build workgraph
wg = build_core_workgraph(
    structures_dir='structures',
    bulk_name='ag2o.cif',
    code_label='VASP-VTST-6.4.3@bohr',
    potential_family='PBE',
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    kpoints_spacing=0.4,
    bulk_parameters={
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NELM': 120,
    },
    bulk_options={
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 40},
        'max_wallclock_seconds': 3600,
        'queue_name': 'normal',
    },
    clean_workdir=False,

    # Metal and oxygen references
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},
    metal_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6},
    oxygen_parameters={'PREC': 'Accurate', 'ENCUT': 520, 'EDIFF': 1e-6},
    metal_options={'resources': {'num_machines': 1}, 'max_wallclock_seconds': 3600},
    oxygen_options={'resources': {'num_machines': 1}, 'max_wallclock_seconds': 3600},

    # Slab generation (will be used as input_slabs in this example)
    miller_indices=[1, 0, 0],
    min_slab_thickness=18.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
    slab_parameters={
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'NSW': 100,
        'IBRION': 2,
    },

    # Bulk electronic properties
    compute_electronic_properties_bulk=True,
    bands_parameters=bulk_defaults,
    bands_options={'resources': {'num_machines': 2}, 'max_wallclock_seconds': 7200},
    band_settings=bulk_defaults['band_settings'],

    # NEW: Slab electronic properties
    compute_electronic_properties_slabs=True,
    slab_electronic_properties={
        'term_0': {
            'bands_parameters': slab_defaults,
            'bands_options': {'resources': {'num_machines': 2}},
            'band_settings': slab_defaults['band_settings'],
        },
        'term_2': {
            'bands_parameters': slab_custom,  # Override with custom
            'bands_options': {'resources': {'num_machines': 4}},  # More resources
            'band_settings': slab_custom['band_settings'],
        },
    },
    slab_bands_parameters=slab_defaults,  # Global defaults
    slab_bands_options={'resources': {'num_machines': 2}},
    slab_band_settings=slab_defaults['band_settings'],

    name='Ag2O_Slab_Electronic_Properties',
)

# Note: For this example to work, you need to first generate slabs, then use them as input_slabs
# This is a conceptual example showing the interface

print(f"WorkGraph built: {wg.name}")
print(f"Number of tasks: {len(wg.tasks)}")
print("\nTo submit:")
print("  from aiida.engine import submit")
print("  node = submit(wg)")
print("  # Wait for completion")
print("  print(node.outputs.slab_bands.keys())  # ['term_0', 'term_2']")
print("  print(node.outputs.slab_dos.term_0)")
