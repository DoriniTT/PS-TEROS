# Surface Hydroxylation Module

Automated workflow for studying surface hydroxylation and vacancy formation on oxide surfaces.

## Overview

This module provides tools for:
- Generating hydroxylated and vacancy-containing surface structures
- Running DFT relaxations on modified surfaces
- Calculating surface energies as a function of temperature and chemical potential
- Analyzing thermodynamic stability of different surface terminations

## Features

- **Automatic Structure Generation**: Generate OH-adsorbed and oxygen-vacancy surfaces from a pristine slab
- **Batch Relaxation**: Efficiently relax multiple surface configurations in parallel
- **Surface Energy Calculations**: Compute surface energies using three reaction pathways
- **JANAF Thermodynamics**: Temperature-dependent thermodynamic data integration
- **Restart Capability**: Resume incomplete calculations from previous runs

## Usage

### Basic Workflow

```python
from aiida import orm
from teros.core.surface_hydroxylation import build_surface_hydroxylation_workgraph

# Define surface modification parameters
surface_params = {
    'oh_coverage': [0.25, 0.5, 0.75, 1.0],  # OH coverage fractions
    'vacancy_coverage': [0.25, 0.5],         # O vacancy fractions
    'sites': 'all',                          # Use all surface sites
}

# Build the workflow
wg = build_surface_hydroxylation_workgraph(
    structure=pristine_slab_structure,
    surface_params=surface_params,
    code=vasp_code,
    builder_inputs={
        'parameters': {'incar': {...}},
        'options': {'resources': {...}},
        'potential_family': 'PBE',
        'potential_mapping': {...},
        'kpoints_spacing': 0.03,
    },
    bulk_structure=bulk_structure,
    bulk_builder_inputs={...},
    max_concurrent_jobs=4,
    calculate_surface_energies=True,
)

wg.submit(wait=False)
```

### With Atom Fixing (Selective Dynamics)

```python
wg = build_surface_hydroxylation_workgraph(
    structure=pristine_slab,
    surface_params=surface_params,
    code=vasp_code,
    builder_inputs=builder_inputs,
    bulk_structure=bulk_structure,
    bulk_builder_inputs=bulk_builder_inputs,

    # Fix bottom layers during relaxation
    fix_type='bottom',
    fix_thickness=5.0,  # Angstroms
    fix_elements=['La', 'Mn'],  # Optional: fix only specific elements
)
```

### Restart from Previous Calculation

```python
wg = build_surface_hydroxylation_workgraph(
    structure=pristine_slab,
    surface_params=surface_params,
    code=vasp_code,
    builder_inputs=builder_inputs,
    bulk_structure=bulk_structure,
    bulk_builder_inputs=bulk_builder_inputs,
    restart_from_pk=12345,  # PK of previous workflow
)
```

## Module Structure

```
surface_hydroxylation/
├── __init__.py                 # Public exports
├── workgraph.py                # Main WorkGraph builder
├── relaxations.py              # Batch relaxation with semaphore control
├── surface_modes.py            # SurfaceModifier class for structure generation
├── surface_energy.py           # Surface energy calculation functions
├── surface_energy_workgraph.py # Surface energy task wrappers
├── thermodynamics.py           # JANAF database integration
├── tasks.py                    # Helper calcfunctions
├── utils.py                    # Utility functions
├── scripts/                    # Data extraction scripts
│   └── extract_janaf_data.py
└── examples/                   # Usage examples
    └── thermodynamics_usage.py
```

## Key Classes and Functions

### `SurfaceModifier`
Generates modified surface structures from a pristine slab:
```python
from teros.core.surface_hydroxylation import SurfaceModifier

modifier = SurfaceModifier(pristine_structure)
oh_structures = modifier.add_oh_groups(coverage=0.5)
vacancy_structures = modifier.create_vacancies(coverage=0.25)
```

### `JanafDatabase`
Access thermodynamic data for temperature-dependent calculations:
```python
from teros.core.surface_hydroxylation import JanafDatabase

db = JanafDatabase()
h2o_data = db.get_species('H2O')
delta_g = db.get_gibbs_energy('H2O', temperature=300)
```

### Surface Energy Functions

Three reaction pathways are supported:

1. **Reaction 1**: H2O adsorption/desorption
   ```
   Surface + H2O ↔ Surface-OH + ½H2
   ```

2. **Reaction 2**: Oxygen vacancy formation
   ```
   Surface-O ↔ Surface-Vo + ½O2
   ```

3. **Reaction 3**: Combined hydroxylation
   ```
   Surface + ½H2O ↔ Surface-OH
   ```

```python
from teros.core.surface_hydroxylation import (
    calc_delta_g_reaction1,
    calc_delta_g_reaction2,
    calc_delta_g_reaction3,
    calculate_surface_energies,
)

# Calculate surface energies for all structures
results = calculate_surface_energies(
    structures=relaxed_structures,
    energies=dft_energies,
    bulk_energy=bulk_energy,
    temperature=300,
    pressure=1.0,
)
```

## Outputs

After workflow completion, access results via:

```python
results = wg.outputs

# Manifest of generated structures
manifest = results.manifest.get_dict()

# Relaxed structures and energies
structures = results.structures  # Dict of StructureData
energies = results.energies      # Dict of Float

# Bulk and pristine references
bulk_structure = results.bulk_structure
bulk_energy = results.bulk_energy
pristine_structure = results.pristine_structure
pristine_energy = results.pristine_energy

# Surface energies for each reaction pathway
gamma_r1 = results.surface_energies_reaction1.get_dict()
gamma_r2 = results.surface_energies_reaction2.get_dict()
gamma_r3 = results.surface_energies_reaction3.get_dict()
```

## Surface Energy Output Format

Each `surface_energies_reactionX` Dict contains:
```python
{
    'temperatures': [300, 400, 500, ...],  # K
    'pressures': [1e-10, 1e-5, 1.0, ...],  # atm
    'gamma': {
        'structure_label': [[...], [...], ...],  # 2D array [T, P]
    },
    'delta_g': {
        'structure_label': [[...], [...], ...],  # 2D array [T, P]
    },
}
```

## Dependencies

- AiiDA-VASP for DFT calculations
- JANAF thermodynamic database (included)
- pymatgen for structure manipulation

## Examples

See `examples/thermodynamics_usage.py` for a complete working example demonstrating:
- Structure generation with different coverages
- Running the full workflow
- Extracting and plotting surface energy diagrams
