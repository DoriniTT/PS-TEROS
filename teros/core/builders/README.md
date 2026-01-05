# Builders Module

Default parameter builders for PS-TEROS workflows.

## Overview

This module provides pre-configured parameter sets for common material systems and calculation types. Instead of manually specifying dozens of VASP parameters, you can use a builder to get sensible defaults and override only what you need.

## Available Builders

| Builder | Description |
|---------|-------------|
| `get_ag2o_defaults` | Binary oxide Ag2O system |
| `get_ag3po4_defaults` | Ternary oxide Ag3PO4 system |
| `get_aimd_defaults` | AIMD parameters for VASP |
| `get_aimd_defaults_cp2k` | AIMD parameters for CP2K |
| `get_electronic_properties_defaults` | DOS/band structure for bulk |
| `get_slab_electronic_properties_defaults` | DOS/band structure for slabs |

## Usage

### Material-Specific Defaults

```python
from teros.core.builders import get_ag2o_defaults
from teros.core import build_core_workgraph

# Get all defaults for Ag2O
defaults = get_ag2o_defaults(
    structures_dir='/path/to/structures',
    code_label='VASP-6.4.1@cluster',
    potential_family='PBE.54',
)

# Override specific parameters
defaults['miller_indices'] = [1, 1, 0]  # (110) surface
defaults['bulk_parameters']['incar']['ENCUT'] = 600

# Build workflow with defaults
wg = build_core_workgraph(**defaults, name='Ag2O_110')
```

### AIMD Parameters (VASP)

```python
from teros.core.builders import get_aimd_defaults

# Get base AIMD parameters (returns INCAR dict)
aimd_params = get_aimd_defaults(
    energy_cutoff=400,              # ENCUT in eV
    electronic_convergence=1e-5,    # EDIFF
    timestep=1.0,                   # POTIM in fs
    smass=0.0,                      # Nosé mass (0=automatic)
    mdalgo=2,                       # NVT Nosé-Hoover thermostat
)

# Use in workflow - temperature and steps are set per stage
wg = build_core_workgraph(
    run_aimd=True,
    aimd_parameters=aimd_params,
    aimd_sequence=[
        {'TEBEG': 300, 'NSW': 500},   # Stage 1: 300K, 500 steps
        {'TEBEG': 600, 'NSW': 1000},  # Stage 2: 600K, 1000 steps
    ],
    ...
)
```

**Note:** Temperature (`TEBEG`) and number of steps (`NSW`) are set in `aimd_sequence`, not in `get_aimd_defaults()`.

### AIMD Parameters (CP2K)

```python
from teros.core.builders import (
    get_aimd_defaults_cp2k,
    get_basis_molopt_content,
    get_gth_potentials_content,
)
from aiida import orm

# Get CP2K AIMD parameters
aimd_params = get_aimd_defaults_cp2k(
    cutoff=400,         # CUTOFF in Ry
    rel_cutoff=60,      # REL_CUTOFF in Ry
    timestep=0.5,       # fs
    eps_scf=1e-6,       # SCF convergence
    thermostat='NOSE',  # NOSE, CSVR, etc.
)

# Create basis and pseudo files (required for CP2K)
basis_file = orm.SinglefileData.from_string(
    get_basis_molopt_content(),
    filename='BASIS_MOLOPT'
)
pseudo_file = orm.SinglefileData.from_string(
    get_gth_potentials_content(),
    filename='GTH_POTENTIALS'
)
```

**Note:** Temperature and steps are set in `aimd_sequence` with keys `'temperature'` and `'steps'` (not VASP-style `TEBEG`/`NSW`).

### Electronic Properties

```python
from teros.core.builders import (
    get_electronic_properties_defaults,
    get_slab_electronic_properties_defaults,
)

# For bulk DOS and band structure calculation
bulk_elec = get_electronic_properties_defaults(
    energy_cutoff=500,              # ENCUT in eV
    kpoints_mesh_density=0.3,       # K-point density
    band_kpoints_distance=0.2,      # K-point spacing for bands
    dos_kpoints_distance=0.2,       # K-point spacing for DOS
    nedos=2000,                     # Number of DOS grid points
    ispin=2,                        # Spin-polarized
)

# For slab DOS (non-periodic in z)
slab_elec = get_slab_electronic_properties_defaults(
    energy_cutoff=500,
    kpoints_mesh_density=0.3,
    nedos=2000,
)
```

**Note:** These return full configuration dicts for the `vasp.v2.bands` workchain including `band_settings` and stage-specific INCAR parameters.

## Builder Structure

Each material builder (e.g., `get_ag2o_defaults`) returns a complete parameter dictionary including:

- `structures_dir` - Path to structure files
- `code_label` - VASP code identifier
- `potential_family` - PAW potential family
- `bulk_parameters` - INCAR for bulk relaxation
- `bulk_options` - Scheduler options for bulk
- `bulk_potential_mapping` - Element to potential mapping
- `slab_parameters` - INCAR for slab relaxation
- `slab_options` - Scheduler options for slabs
- `miller_indices` - Surface orientation
- `min_slab_thickness` - Minimum slab thickness
- `min_vacuum_thickness` - Vacuum gap size
- Reference parameters (metal, nonmetal, oxygen)

## Creating Custom Builders

To create a builder for a new material system:

```python
from copy import deepcopy
from teros.core.builders.default_ag2o_builders import update_builder_params

def get_my_material_defaults(
    structures_dir=None,
    code_label=None,
    potential_family=None,
    **overrides
):
    """Get defaults for my material system."""

    defaults = {
        'structures_dir': structures_dir,
        'code_label': code_label,
        'potential_family': potential_family,

        'bulk_parameters': {
            'incar': {
                'ENCUT': 520,
                'EDIFF': 1e-6,
                # ... more parameters
            }
        },
        # ... more sections
    }

    # Apply any overrides
    if overrides:
        defaults = update_builder_params(defaults, overrides)

    return defaults
```

## Module Structure

- `__init__.py` - Public exports
- `default_ag2o_builders.py` - Ag2O binary oxide defaults
- `default_ag3po4_builders.py` - Ag3PO4 ternary oxide defaults
- `aimd_builder.py` - VASP AIMD parameter builder
- `aimd_builder_cp2k.py` - CP2K AIMD parameter builder
- `electronic_properties_builder.py` - DOS/band structure builders
