# Working with PS-TEROS and AiiDA-WorkGraph

This is a comprehensive guide for developing and working with the **PS-TEROS** code using **AiiDA**, **AiiDA-WorkGraph**, and **AiiDA-VASP**.

---

## Quick Start Commands

### Before you begin

Make sure that you are in the 'presto' profile in AiiDA by doing `verdi profile set-default presto`, then do `verdi status` to check if everything is working fine.
If the daemon is not running, you can start it by doing `verdi daemon start`.

### After every modification

```bash
verdi daemon restart
```

### Python binary

```bash
source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python
```

### Launch workflow

```bash
source ~/envs/aiida/bin/activate && /home/thiagotd/envs/aiida/bin/python /home/thiagotd/git/PS-TEROS/examples/vasp/step_x_example.py
```

### Wait for results

Use `sleep 15` (or more) before analyzing nodes.

### Analyze results

```bash
verdi process show <PK>
verdi process report <PK>
```

## Testing Configuration

### Test Localhost (vasp-6.5.1-std@localhost)

```
Code Label: vasp-6.5.1-std@localhost
PK: 20510
Processors: 8
Queue: None (no queue system)
```

```python
options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 8,
    },
}
# IMPORTANT: Only one job at a time
max_concurrent_jobs = 1
```

### Error Recovery

When fixing errors, kill all processes:

```bash
# Kill AiiDA workflow
verdi process kill <WORKGRAPH_PK>

# Kill remote VASP processes
killall vasp_std

# Clear cache and restart daemon
find . -type d -name __pycache__ -exec rm -rf {} +
verdi daemon restart
```

---

## Validation Checklist

Before merging a new feature:

- [ ] Main node returns exit code `[0]`
- [ ] Feature outputs appear correctly in `verdi process show <PK>`
- [ ] Test with low precision parameters (proper ENCUT, NSW, EDIFF, etc.)
- [ ] No Python cache issues (run cache clear)
- [ ] Example script follows naming convention `step_XX_feature_name.py`
- [ ] Module exported in `__init__.py`
- [ ] Documentation in example folder (if needed)

---

## Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| Import errors after changes | Run `verdi daemon restart` |
| Tasks not running in parallel | Check `max_number_jobs` setting |
| Serialization errors | Use AiiDA types (orm.Float, orm.Dict) not Python types |
| Provenance cycles | Use `@task.graph` instead of `@task.calcfunction` for stored nodes |
| Remote folder not cleaning | Set `clean_workdir=True` in builder inputs |

---

## Branch Management

1. **Development**: Create feature folders in `teros/experimental/`
2. **Testing**: Create test scripts in `examples/`
3. **Integration**: Move validated code to `teros/core/`
4. **Documentation**: Create docs alongside test scripts

---

## Useful Commands Reference

```bash
# Profile management
verdi profile list
verdi profile set-default presto

# Process monitoring
verdi process status <PK>
verdi process show <PK>
verdi process report <PK>
verdi process kill <PK>

# Node inspection
verdi node show <PK>
verdi data core.dict show <PK>
verdi data core.structure show <PK>

# Daemon management
verdi daemon status
verdi daemon start
verdi daemon restart
verdi daemon stop

# Code information
verdi code show <CODE_LABEL>
verdi code list
```

### Clear Python cache

```bash
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
```

---

## Repository Structure

```
/home/thiagotd/git/PS-TEROS/
├── teros/                    # Main package
│   ├── core/                 # Core modules (main development)
│   │   ├── __init__.py       # Public API exports
│   │   ├── workgraph.py      # Main WorkGraph entry point
│   │   ├── workflow_presets.py # Workflow preset system
│   │   ├── helper_functions.py # Generic helper functions
│   │   ├── slabs.py          # Slab generation and relaxation
│   │   ├── thermodynamics.py # Surface energy calculations
│   │   ├── cleavage.py       # Cleavage energy calculations
│   │   ├── hf.py             # Formation enthalpy
│   │   ├── adsorption_energy.py # Adsorption energy module
│   │   ├── fixed_atoms.py    # Atom constraint utilities
│   │   ├── aimd_functions.py # AIMD for VASP
│   │   ├── aimd_cp2k.py      # AIMD for CP2K
│   │   ├── builders/         # Builder utilities
│   │   ├── custom_calculation/ # Custom calculation support
│   │   ├── surface_hydroxylation/ # Surface hydroxylation module
│   │   └── tests/            # Unit tests
│   └── experimental/         # New features under development
├── examples/                 # Test scripts organized by feature
│   ├── vasp/                 # VASP-based examples (step_01 to step_19)
│   ├── cp2k/                 # CP2K-based examples
│   └── structures/           # CIF/POSCAR files
├── docs/                     # Documentation
└── setup.py                  # Package installation
```

---

## Core Module Architecture

### Key Design Patterns

PS-TEROS uses **AiiDA-WorkGraph** for workflow orchestration with these patterns:

1. **Scatter-Gather Pattern**: Process multiple structures in parallel
2. **Workflow Presets**: Pre-configured calculation profiles
3. **Builder Inputs**: Fine-grained control over calculation parameters

### Decorator Types

| Decorator | Use Case | Provenance |
|-----------|----------|------------|
| `@task.calcfunction` | Pure Python functions, small computations | Creates AiiDA CalcFunctionNode |
| `@task.graph` | Subworkflows, parallel execution | Creates nested WorkGraph |
| `@task(WorkChain)` | Wrap existing WorkChains (VASP, CP2K) | Uses existing WorkChain |

### Entry Point: `build_core_workgraph()`

The main entry point in `workgraph.py` builds workflows using presets:

```python
from teros.core.workgraph import build_core_workgraph

wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # or bulk_only, aimd_only, etc.
    structures_dir='structures/',
    bulk_name='ag2o.cif',
    # ... other parameters
)
wg.submit(wait=False)
```

### Available Workflow Presets

| Preset | Description |
|--------|-------------|
| `bulk_only` | Bulk relaxation only |
| `formation_enthalpy_only` | Formation enthalpy without surfaces |
| `surface_thermodynamics` | Complete surface energy workflow (default) |
| `cleavage_only` | Cleavage energy calculations |
| `relaxation_energy_only` | Relaxation energy calculations |
| `electronic_structure_bulk_only` | DOS/bands for bulk |
| `electronic_structure_slabs_only` | DOS/bands for slabs |
| `aimd_only` | AIMD simulation |
| `adsorption_energy` | Adsorption energy calculations |
| `comprehensive` | All features enabled |

---

## Creating New Modules

### Step 1: Create the Module File

Create a new file in `teros/core/` following this template:

```python
"""
Module Description

This module provides [description of functionality].
"""

from __future__ import annotations

import typing as t
from aiida import orm
from aiida_workgraph import task, namespace, dynamic


# =============================================================================
# CALCFUNCTIONS (Pure Python, creates provenance)
# =============================================================================

@task.calcfunction
def my_calculation(
    structure: orm.StructureData,
    energy: orm.Float,
    parameters: orm.Dict,
) -> orm.Dict:
    """
    Perform calculation on a single structure.

    All inputs MUST be AiiDA data types (orm.StructureData, orm.Float, etc.)
    Return type MUST be AiiDA data type or dict annotation.

    Args:
        structure: Input structure
        energy: Energy value
        parameters: Calculation parameters

    Returns:
        Dictionary with results
    """
    # Access values with .value or .get_dict()
    struct_ase = structure.get_ase()
    energy_val = energy.value
    params = parameters.get_dict()

    # Perform calculation
    result = {...}

    # Return AiiDA Dict
    return orm.Dict(dict=result)


# =============================================================================
# TASK GRAPHS (Parallel execution, scatter-gather)
# =============================================================================

@task.graph
def my_scatter_function(
    structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    code: orm.Code,
    parameters: t.Mapping[str, t.Any],
    options: t.Mapping[str, t.Any],
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(
    results=dynamic(orm.Dict),
    energies=dynamic(orm.Float),
)]:
    """
    Scatter-gather pattern for parallel processing.

    Args:
        structures: Dynamic namespace of structures to process
        code: AiiDA code for calculations
        parameters: Calculation parameters
        options: Scheduler options
        max_number_jobs: Concurrency limit

    Returns:
        Dictionary with output namespaces
    """
    from aiida.plugins import WorkflowFactory
    from aiida_workgraph import get_current_graph

    # Set concurrency limit
    if max_number_jobs is not None:
        wg = get_current_graph()
        wg.max_number_jobs = max_number_jobs

    # Load WorkChain
    VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
    vasp_task = task(VaspWorkChain)

    # Output dictionaries
    results_dict: dict[str, orm.Dict] = {}
    energies_dict: dict[str, orm.Float] = {}

    # Scatter: process each structure in parallel
    for label, structure in structures.items():
        # Build inputs
        calc_inputs = {
            'structure': structure,
            'code': code,
            'parameters': {'incar': dict(parameters)},
            'options': dict(options),
        }

        # Create task
        calc = vasp_task(**calc_inputs)

        # Collect outputs
        results_dict[label] = calc.misc
        energies_dict[label] = extract_energy(calc.misc).result

    # Gather: return collected results
    return {
        'results': results_dict,
        'energies': energies_dict,
    }
```

### Step 2: Export in `__init__.py`

Add your functions to `teros/core/__init__.py`:

```python
from .my_module import (
    my_calculation,
    my_scatter_function,
)

__all__ = [
    # ... existing exports
    'my_calculation',
    'my_scatter_function',
]
```

### Step 3: Integrate with Workgraph (Optional)

If your module should be part of the main workflow, modify `workgraph.py`:

1. Add parameters to `build_core_workgraph()` signature
2. Add corresponding flags to `workflow_presets.py`
3. Wire the tasks in `core_workgraph()` or `build_core_workgraph()`

### Step 4: Create Test Script

Create a test in `examples/vasp/step_XX_my_feature.py`:

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
STEP XX: My Feature Test

Description of what this tests.
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

def main():
    print("\n" + "="*70)
    print("STEP XX: MY FEATURE")
    print("="*70)

    load_profile(profile='psteros')

    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    wg = build_core_workgraph(
        workflow_preset='my_preset',
        structures_dir=structures_dir,
        # ... parameters
        name='StepXX_MyFeature',
    )

    wg.submit(wait=False)

    print(f"WorkGraph PK: {wg.pk}")
    print(f"Monitor: verdi process status {wg.pk}")

    return wg

if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
```

---

## Important AiiDA-WorkGraph Patterns

### Dynamic Namespaces

Use for dict outputs where keys are determined at runtime:

```python
# Input annotation
structures: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)]

# Output annotation
-> t.Annotated[dict, namespace(results=dynamic(orm.Dict))]
```

### Task Dependencies

Chain tasks sequentially:

```python
task_a >> task_b  # task_b waits for task_a
group(task_a, task_b) >> task_c  # task_c waits for both
```

### Accessing Task Outputs

```python
# Socket reference (for wiring)
output_socket = task.outputs.structure

# Result reference (for using in graph)
result = my_calcfunction(data=other_task.misc).result
```

### Deep Merge for Parameter Overrides

```python
from copy import deepcopy

def deep_merge_dicts(base: dict, override: dict) -> dict:
    """Recursively merge override into base."""
    result = deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge_dicts(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result
```

---

