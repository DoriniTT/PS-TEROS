# Custom Calculation Module

Generic VASP calculation wrapper for PS-TEROS workflows.

## Overview

This module provides a flexible way to run arbitrary VASP calculations outside of the predefined workflow presets. It supports both single and multiple structure calculations with optional selective dynamics (atom fixing).

## Features

- Single structure or batch calculations
- Automatic energy extraction
- Optional selective dynamics (fix atoms by position or element)
- Concurrency control for batch calculations
- Direct access to VASP outputs (structure, energy, misc)

## Usage

### Single Structure Calculation

```python
from aiida import orm
from teros.core.custom_calculation import build_custom_calculation_workgraph

# Load your structure
structure = orm.StructureData(ase=your_ase_structure)

# Define VASP parameters
builder_inputs = {
    'parameters': {
        'incar': {
            'ENCUT': 520,
            'EDIFF': 1e-6,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'NSW': 100,
            'ISIF': 2,
        }
    },
    'options': {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
        'queue_name': 'normal',
        'walltime_seconds': 86400,
    },
    'kpoints_spacing': 0.03,
    'potential_family': 'PBE.54',
    'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
    'clean_workdir': False,
}

# Build and submit
wg = build_custom_calculation_workgraph(
    structure=structure,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,
    name='my_calculation',
)
wg.submit(wait=False)
```

### Multiple Structures

```python
structures = [structure1, structure2, structure3]

wg = build_custom_calculation_workgraph(
    structure=structures,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,  # Same params for all
    max_concurrent_jobs=2,          # Limit parallelism
)
```

### With Selective Dynamics

```python
wg = build_custom_calculation_workgraph(
    structure=slab_structure,
    code_label='VASP-6.4.1@cluster',
    builder_inputs=builder_inputs,
    fix_type='bottom',        # Fix bottom atoms
    fix_thickness=3.0,        # 3 Angstroms from bottom
    fix_elements=['Ag'],      # Only fix Ag atoms (optional)
)
```

### Extracting Results

```python
from teros.core.custom_calculation import get_custom_results

# After workflow completes
results = get_custom_results(wg)

print(f"Energy: {results['energies']} eV")
print(f"Structure: {results['structures']}")
print(f"VASP misc: {results['misc']}")
```

## Module Structure

- `__init__.py` - Public exports
- `workgraph.py` - WorkGraph builder and result extraction
- `tasks.py` - Helper calcfunctions (energy extraction)

## API Reference

### `build_custom_calculation_workgraph`

Main builder function for custom VASP calculations.

**Parameters:**
- `structure`: Single `StructureData` or list of structures
- `code_label`: VASP code label (e.g., `'VASP-6.4.1@cluster'`)
- `builder_inputs`: Dict with VASP parameters (see example above)
- `name`: WorkGraph name (default: `'custom_calc'`)
- `max_concurrent_jobs`: Limit concurrent calculations (optional)
- `fix_type`: `'bottom'`, `'top'`, `'center'`, or `None`
- `fix_thickness`: Thickness in Angstroms for fixing region
- `fix_elements`: List of element symbols to fix (optional)

**Returns:** WorkGraph ready to submit

### `get_custom_results`

Extract results from completed WorkGraph.

**Parameters:**
- `workgraph`: Completed WorkGraph

**Returns:** Dict with `energies`, `structures`, and `misc` keys
