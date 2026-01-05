# PS-TEROS Developer Guide

A comprehensive guide for developing and working with **PS-TEROS** (Predicting Stability of TERminations of Oxide Surfaces) using **AiiDA** and **AiiDA-WorkGraph**.
This is a comprehensive guide for developing and working with the **PS-TEROS** code using **AiiDA**, **AiiDA-WorkGraph**, and **AiiDA-VASP**.

## Table of Contents

1. [Overview](#1-overview)
2. [Environment Setup](#2-environment-setup)
3. [Project Structure](#3-project-structure)
4. [Core Architecture Concepts](#4-core-architecture-concepts)
5. [Building New Modules](#5-building-new-modules)
6. [Updating Existing Modules](#6-updating-existing-modules)
7. [VASP + AiiDA-VASP Integration](#7-vasp--aiida-vasp-integration)
8. [AiiDA Commands Reference](#8-aiida-commands-reference)
9. [Testing & Validation](#9-testing--validation)
10. [Troubleshooting](#10-troubleshooting)
11. [Code Style & Conventions](#11-code-style--conventions)

---

## 1. Overview

### What is PS-TEROS?

PS-TEROS is a Python package for automated computational materials science workflows, specializing in:

- **Surface thermodynamics** - Surface energy calculations for oxides and metals
- **Formation enthalpy** - Standard enthalpies of formation
- **Cleavage energy** - Energy to split crystals into complementary surfaces
- **Electronic properties** - DOS and band structure calculations
- **AIMD** - Ab initio molecular dynamics simulations
- **Adsorption energy** - Molecule-surface interaction calculations
- **Surface hydroxylation** - OH/vacancy generation and thermodynamics

### Technology Stack

| Component | Purpose |
|-----------|---------|
| **AiiDA** | Workflow management, provenance tracking, data persistence |
| **AiiDA-WorkGraph** | Task-based workflow composition with scatter-gather patterns |
| **AiiDA-VASP** | VASP calculation interface for DFT calculations |
| **Pymatgen** | Structure manipulation, slab generation |
| **ASE** | Atomic structure I/O |

---

## 2. Environment Setup

### Before You Begin

1. **Activate the AiiDA environment:**
   ```bash
   source ~/envs/aiida/bin/activate
   ```

2. **Set the correct AiiDA profile:**
   ```bash
   verdi profile set-default presto
   ```

3. **Verify AiiDA status:**
   ```bash
   verdi status
   ```

   Expected output shows all services running:
   ```
   ✔ version:     AiiDA v2.x.x
   ✔ config:      /home/user/.aiida
   ✔ profile:     presto
   ✔ storage:     SqliteDosStorage[...]
   ✔ rabbitmq:    Connected
   ✔ daemon:      Running
   ```

4. **Start the daemon if not running:**
   ```bash
   verdi daemon start
   ```

### After Every Code Modification

**Always restart the daemon** to reload Python modules:
```bash
verdi daemon restart
```

### Python Binary & Scripts

```bash
# Python interpreter
source ~/envs/aiida/bin/activate && python

# Run a workflow script
source ~/envs/aiida/bin/activate && python /path/to/script.py
```

### Clear Python Cache (when needed)

```bash
find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete
```

---

## 3. Project Structure

```
PS-TEROS/
├── teros/                              # Main package
│   ├── core/                           # Production code
│   │   ├── __init__.py                 # Module exports (public API)
│   │   ├── workgraph.py                # Main orchestration workflow (~2000 lines)
│   │   ├── workflow_presets.py         # Preset system (11 presets)
│   │   │
│   │   ├── # Physics Modules
│   │   ├── thermodynamics.py           # Surface energy γ(Δμ_O, Δμ_M)
│   │   ├── hf.py                       # Formation enthalpy
│   │   ├── cleavage.py                 # Cleavage energy (complementary pairs)
│   │   ├── slabs.py                    # Slab generation & relaxation
│   │   ├── adsorption_energy.py        # Adsorption calculations
│   │   ├── fixed_atoms.py              # Selective dynamics (atom fixing)
│   │   │
│   │   ├── # Specialized Submodules
│   │   ├── aimd/                       # Ab initio molecular dynamics
│   │   ├── custom_calculation/         # Generic VASP calculations
│   │   ├── metal_surface_energy/       # Elemental metal surfaces
│   │   ├── surface_hydroxylation/      # OH/vacancy studies
│   │   │
│   │   ├── builders/                   # Default parameter builders
│   │   │   ├── default_ag2o_builders.py
│   │   │   ├── default_ag3po4_builders.py
│   │   │   ├── aimd_builder.py
│   │   │   └── electronic_properties_builder.py
│   │   │
│   │   └── tests/                      # Unit tests
│   │
│   └── experimental/                   # New features under development
│
├── examples/                           # Test scripts & usage examples
│   ├── vasp/                           # Main example series (step_01 → step_18)
│   │   ├── step_01_bulk_only.py
│   │   ├── step_02_formation_enthalpy.py
│   │   ├── ...
│   │   └── structures/                 # CIF files
│   ├── custom_calculation/
│   ├── surface_thermo_gold/
│   └── ...
│
├── docs/                               # Documentation
├── legacy/                             # Archived code
├── CLAUDE.md                           # This file
├── CHANGE.md                           # Changelog
└── setup.py                            # Installation
```

### Module Organization Pattern

Each specialized submodule follows this structure:
```
module_name/
├── __init__.py         # Public exports
├── workgraph.py        # Main WorkGraph builder function
├── tasks.py            # @task.calcfunction helpers
├── utils.py            # Pure Python utilities (optional)
└── README.md           # Module documentation
```

---

## 4. Core Architecture Concepts

### 4.1 Task Decorators

All executable functions use `@task` decorators from `aiida_workgraph`:

```python
from aiida import orm
from aiida_workgraph import task

# Calcfunction - single computation node (stored in database)
@task.calcfunction
def extract_total_energy(misc: orm.Dict) -> orm.Float:
    """Pure computation, inputs/outputs tracked by AiiDA."""
    energy = misc.get_dict()['total_energies']['energy_extrapolated']
    return orm.Float(energy)

# Graph task - composes multiple tasks into a sub-workflow
@task.graph(outputs=['energy', 'structure'])
def my_workflow(structure: orm.StructureData, ...) -> dict:
    """Orchestrates multiple tasks."""
    vasp_task = VaspTask(structure=structure, ...)
    energy = extract_total_energy(misc=vasp_task.outputs.misc)
    return {'energy': energy, 'structure': vasp_task.outputs.structure}
```

### 4.2 Wrapping External WorkChains

```python
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task

# Load the VASP WorkChain
VaspWorkChain = WorkflowFactory('vasp.v2.vasp')

# Wrap it as a task for use in WorkGraphs
VaspTask = task(VaspWorkChain)

# Now use in a workflow
vasp_task = VaspTask(
    structure=structure,
    code=code,
    parameters=parameters,
    ...
)
```

### 4.3 Dynamic Namespaces (Scatter-Gather Pattern)

For parallel execution of variable-length collections:

```python
from aiida_workgraph import dynamic, namespace
import typing as t

@task.graph(outputs=['relaxed_structures', 'energies'])
def relax_slabs_scatter(
    # Input: dynamic dict of structures
    slab_structures: t.Annotated[dict, dynamic(orm.StructureData)],
    ...
) -> t.Annotated[dict, namespace(
    relaxed_structures=dynamic(orm.StructureData),
    energies=dynamic(orm.Float),
)]:
    """
    Relax each slab in parallel.

    Input: {'term_0': StructureData, 'term_1': StructureData, ...}
    Output: {'term_0': Float, 'term_1': Float, ...}
    """
    outputs = {'relaxed_structures': {}, 'energies': {}}

    for label, slab in slab_structures.items():
        # Each creates a separate task node
        vasp_task = VaspTask(structure=slab, name=f'relax_{label}', ...)
        energy = extract_total_energy(misc=vasp_task.outputs.misc)

        outputs['relaxed_structures'][label] = vasp_task.outputs.structure
        outputs['energies'][label] = energy.outputs.result

    return outputs
```

### 4.4 Builder Input Pattern

All VASP calculations use a standardized `builder_inputs` dictionary:

```python
builder_inputs = {
    # VASP INCAR parameters
    'parameters': {
        'incar': {
            'ENCUT': 520,
            'EDIFF': 1e-6,
            'ISMEAR': 0,
            'SIGMA': 0.05,
            'IBRION': 2,
            'NSW': 100,
            'ISIF': 3,      # 3=full relax, 2=ions only
            'PREC': 'Accurate',
            'LREAL': 'Auto',
            'ALGO': 'Normal',
            'LWAVE': False,
            'LCHARG': False,
        }
    },

    # Scheduler options
    'options': {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
        'queue_name': 'normal',
        'walltime_seconds': 86400,
        'max_memory_kb': 4000000,
    },

    # K-point density
    'kpoints_spacing': 0.03,  # Å⁻¹

    # Pseudopotentials
    'potential_family': 'PBE.54',
    'potential_mapping': {'Ag': 'Ag', 'O': 'O'},

    # Remote data handling
    'clean_workdir': False,  # Keep files for restart
}
```

### 4.5 Deep Merge for Parameter Overrides

```python
from teros.core.slabs import deep_merge_dicts

base = {'parameters': {'incar': {'ENCUT': 400, 'ISMEAR': 0}}}
override = {'parameters': {'incar': {'ENCUT': 520}}}  # Only change ENCUT

result = deep_merge_dicts(base, override)
# {'parameters': {'incar': {'ENCUT': 520, 'ISMEAR': 0}}}  # ISMEAR preserved!
```

### 4.6 Workflow Presets

High-level convenience with fine-grained control:

```python
from teros.core import build_core_workgraph

# Simple: use a preset
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    structures_dir='/path/to/structures',
    ...
)

# Advanced: override specific flags
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,   # Enable optional feature
    max_concurrent_jobs=4,   # Control parallelism
    ...
)
```

**Available presets:** `surface_thermodynamics`, `bulk_only`, `cleavage_only`, `relaxation_energy_only`, `electronic_structure_bulk_only`, `electronic_structure_slabs_only`, `aimd_only`, `comprehensive`

---

## 5. Building New Modules

### Step 1: Create Module Structure

```bash
mkdir -p teros/core/my_new_module
touch teros/core/my_new_module/__init__.py
touch teros/core/my_new_module/workgraph.py
touch teros/core/my_new_module/tasks.py
touch teros/core/my_new_module/README.md
```

### Step 2: Implement Tasks (`tasks.py`)

```python
"""Helper tasks for my_new_module."""

from aiida import orm
from aiida_workgraph import task


@task.calcfunction
def my_calculation(
    structure: orm.StructureData,
    parameters: orm.Dict,
) -> orm.Dict:
    """
    Perform a calculation.

    Args:
        structure: Input structure
        parameters: Calculation parameters

    Returns:
        Results dictionary
    """
    # Get Python objects from AiiDA nodes
    params = parameters.get_dict()

    # Do computation (NO I/O here - pure function)
    result = {'computed_value': 42.0}

    # Return AiiDA node
    return orm.Dict(result)
```

### Step 3: Implement WorkGraph Builder (`workgraph.py`)

```python
"""WorkGraph builder for my_new_module."""

import typing as t
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task, WorkGraph, dynamic, namespace

from .tasks import my_calculation

# Wrap VASP WorkChain
VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
VaspTask = task(VaspWorkChain)


@task.graph(outputs=['result', 'energy'])
def my_module_workgraph(
    structure: orm.StructureData,
    code: orm.Code,
    parameters: orm.Dict,
    options: orm.Dict,
    potential_family: str,
    potential_mapping: dict,
    kpoints_spacing: float = 0.03,
    clean_workdir: bool = False,
) -> dict:
    """
    Main workflow for my_new_module.

    Args:
        structure: Input structure
        code: VASP code
        parameters: VASP INCAR parameters
        options: Scheduler options
        potential_family: PAW potential family
        potential_mapping: Element to potential mapping
        kpoints_spacing: K-point density
        clean_workdir: Whether to clean remote directory

    Returns:
        Dictionary with 'result' and 'energy' outputs
    """
    # Task 1: Run VASP calculation
    vasp_task = VaspTask(
        structure=structure,
        code=code,
        parameters=parameters,
        options=options,
        potential_family=orm.Str(potential_family),
        potential_mapping=orm.Dict(potential_mapping),
        kpoints=orm.KpointsData().set_cell_from_structure(structure).set_kpoints_mesh_from_density(kpoints_spacing),
        clean_workdir=orm.Bool(clean_workdir),
    )

    # Task 2: Extract and process results
    result = my_calculation(
        structure=vasp_task.outputs.structure,
        parameters=vasp_task.outputs.misc,
    )

    return {
        'result': result.outputs.result,
        'energy': vasp_task.outputs.misc,
    }


def build_my_module_workgraph(
    structure_path: str,
    code_label: str,
    builder_inputs: dict,
    name: str = 'MyModule',
    max_concurrent_jobs: int = None,
) -> WorkGraph:
    """
    Build a WorkGraph for my_new_module.

    This is the main user-facing function.

    Args:
        structure_path: Path to CIF/POSCAR file
        code_label: AiiDA code label (e.g., 'VASP-6.5.1@localwork')
        builder_inputs: VASP configuration dictionary
        name: Workflow name
        max_concurrent_jobs: Maximum concurrent VASP calculations

    Returns:
        WorkGraph ready to submit

    Example:
        >>> wg = build_my_module_workgraph(
        ...     structure_path='/path/to/structure.cif',
        ...     code_label='VASP-6.5.1@localwork',
        ...     builder_inputs={...},
        ... )
        >>> wg.submit(wait=False)
    """
    from ase.io import read

    # Load structure
    ase_structure = read(structure_path)
    structure = orm.StructureData(ase=ase_structure)

    # Get code
    code = orm.load_code(code_label)

    # Prepare AiiDA nodes from builder_inputs
    parameters = orm.Dict(builder_inputs.get('parameters', {}))
    options = orm.Dict(builder_inputs.get('options', {}))

    # Build workgraph
    wg = my_module_workgraph(
        structure=structure,
        code=code,
        parameters=parameters,
        options=options,
        potential_family=builder_inputs.get('potential_family', 'PBE'),
        potential_mapping=builder_inputs.get('potential_mapping', {}),
        kpoints_spacing=builder_inputs.get('kpoints_spacing', 0.03),
        clean_workdir=builder_inputs.get('clean_workdir', False),
    )

    wg.name = name

    # Set concurrency limit
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg
```

### Step 4: Export in `__init__.py`

```python
"""My New Module - Brief description."""

from .tasks import my_calculation
from .workgraph import my_module_workgraph, build_my_module_workgraph

__all__ = [
    'my_calculation',
    'my_module_workgraph',
    'build_my_module_workgraph',
]
```

### Step 5: Add to Core Exports

Edit `teros/core/__init__.py`:
```python
from .my_new_module import (
    my_calculation,
    build_my_module_workgraph,
)

__all__ = [
    # ... existing exports ...
    'my_calculation',
    'build_my_module_workgraph',
]
```

### Step 6: Create Example Script

Create `examples/my_new_module/test_my_module.py`:
```python
#!/usr/bin/env python
"""Test script for my_new_module."""

from aiida import load_profile
from teros.core.my_new_module import build_my_module_workgraph

def main():
    load_profile('presto')

    builder_inputs = {
        'parameters': {'incar': {'ENCUT': 520, 'EDIFF': 1e-6, ...}},
        'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        'potential_family': 'PBE',
        'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
        'kpoints_spacing': 0.03,
    }

    wg = build_my_module_workgraph(
        structure_path='structures/ag2o.cif',
        code_label='VASP-6.5.1@localwork',
        builder_inputs=builder_inputs,
        max_concurrent_jobs=1,  # IMPORTANT for localwork
    )

    wg.submit(wait=False)
    print(f"Submitted: PK={wg.pk}")
    print(f"Monitor with: verdi process show {wg.pk}")

if __name__ == '__main__':
    main()
```

---

## 6. Updating Existing Modules

### Adding a New Parameter

1. **Add to function signature:**
   ```python
   def build_core_workgraph(
       ...,
       my_new_param: bool = False,  # NEW: Add with default
   ):
   ```

2. **Use conditionally:**
   ```python
   if my_new_param:
       new_task = my_new_function(...)
   ```

3. **Document in docstring:**
   ```python
   """
   Args:
       my_new_param: Enable new feature (default: False)
   """
   ```

4. **Update CHANGE.md**

### Adding a New Output

1. **Add to `@task.graph(outputs=[...])`:**
   ```python
   @task.graph(outputs=[..., 'new_output'])
   def core_workgraph(...):
   ```

2. **Return in the output dict:**
   ```python
   return {
       ...,
       'new_output': new_task.outputs.result,
   }
   ```

### Modifying VASP Parameters

- Always use `deep_merge_dicts` for parameter overrides
- Never mutate the original dictionary
- Validate critical parameters (e.g., ISIF=3 for bulk)

---

## 7. VASP + AiiDA-VASP Integration

### How VASP Calculations Work in AiiDA

```
User Script → build_workgraph() → WorkGraph.submit()
                                        ↓
                                   AiiDA Daemon
                                        ↓
                              VaspWorkChain (vasp.v2.vasp)
                                        ↓
                              ┌─────────────────────┐
                              │   Remote Cluster    │
                              │                     │
                              │   VASP Execution    │
                              │   ├── INCAR         │
                              │   ├── POSCAR        │
                              │   ├── POTCAR        │
                              │   ├── KPOINTS       │
                              │   └── vasp_std      │
                              │                     │
                              └─────────────────────┘
                                        ↓
                              Parse OUTCAR, vasprun.xml
                                        ↓
                              Store results in AiiDA
```

### VaspWorkChain Outputs

| Output | Type | Description |
|--------|------|-------------|
| `structure` | `StructureData` | Relaxed structure (if NSW > 0) |
| `misc` | `Dict` | Parsed VASP results (energies, forces, etc.) |
| `remote_folder` | `RemoteData` | Link to files on cluster |
| `retrieved` | `FolderData` | Downloaded output files |

### Accessing Energy from VASP

```python
# misc contains parsed VASP output
misc_dict = vasp_task.outputs.misc.get_dict()

# Energy is in total_energies
energy = misc_dict['total_energies']['energy_extrapolated']  # Recommended
# or
energy = misc_dict['total_energies']['energy_no_entropy']
# or
energy = misc_dict['total_energies']['energy']
```

### Important INCAR Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `ENCUT` | Energy cutoff (eV) | 400-600 |
| `EDIFF` | Electronic convergence (eV) | 1e-5 to 1e-7 |
| `EDIFFG` | Ionic convergence (eV/Å) | -0.01 to -0.05 |
| `ISMEAR` | Smearing method | 0 (Gaussian), 1 (Methfessel-Paxton) |
| `SIGMA` | Smearing width (eV) | 0.05-0.2 |
| `NSW` | Max ionic steps | 0 (SCF only), 100+ (relaxation) |
| `IBRION` | Ionic relaxation method | 2 (CG), -1 (no relaxation) |
| `ISIF` | Stress tensor control | 2 (ions), 3 (ions+cell) |

### Metals vs. Insulators

```python
# For insulators/semiconductors
parameters = {'ISMEAR': 0, 'SIGMA': 0.05}

# For metals
parameters = {'ISMEAR': 1, 'SIGMA': 0.2}  # Methfessel-Paxton
```

---

## 8. AiiDA Commands Reference

### Workflow Monitoring

```bash
# Show workflow status (most useful!)
verdi process show <PK>

# Show hierarchical report with sub-processes
verdi process report <PK>

# Watch status in real-time
verdi process status <PK>

# List recent processes
verdi process list -a -p 1  # All processes from past 1 day

# List only running processes
verdi process list
```

### Understanding Process States

| State | Meaning |
|-------|---------|
| `Created` | Process created but not yet submitted |
| `Waiting` | Waiting for sub-processes to complete |
| `Running` | Currently executing |
| `Finished [0]` | Completed successfully |
| `Finished [1-n]` | Completed with error (exit code) |
| `Excepted` | Crashed with exception |
| `Killed` | Manually terminated |

### Process Control

```bash
# Kill a process and all children
verdi process kill <PK>

# Pause a process
verdi process pause <PK>

# Resume a paused process
verdi process play <PK>
```

### Inspecting Results

```bash
# Show all outputs of a node
verdi node show <PK>

# Show attributes of a Dict node
verdi data core.dict show <PK>

# Show structure information
verdi data core.structure show <PK>

# Export structure to file
verdi data core.structure export <PK> -o structure.cif
```

### Code & Computer Management

```bash
# List available codes
verdi code list

# Show code details
verdi code show <CODE_LABEL>

# List computers
verdi computer list

# Test computer connection
verdi computer test <COMPUTER_NAME>
```

### Daemon Control

```bash
verdi daemon start
verdi daemon stop
verdi daemon restart   # IMPORTANT: after code changes!
verdi daemon status
```

### Database Queries

```bash
# Find nodes by type
verdi node find -T WORKGRAPH

# Find by label
verdi node find -l "Au_surface*"

# Show node provenance graph
verdi node graph generate <PK>
```

---

## 9. Testing & Validation

### Test Cluster Configuration

For testing, use `VASP-6.5.1@localwork`:

```python
code_label = 'VASP-6.5.1@localwork'

options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 16,
    },
}
builder.metadata.options.max_wallclock_seconds = 3 * 3600  # Lets limit the calculations to 3 hour. This may be fine because we will use this cluster mainly to do testing.

# CRITICAL: the local work computer can only run ONE job at a time
max_concurrent_jobs = 1
```

### Validation Checklist

Your implementation is ready to merge when:

1. **Main node returns `[0]`** (successful completion)
   ```bash
   verdi process show <PK>
   # Look for: "Finished [0]"
   ```

2. **Expected outputs are present:**
   ```bash
   verdi process show <PK>
   # Check "Outputs" section lists your new outputs
   ```

3. **Values are physically reasonable:**
   - Surface energies: 0.5-3.0 J/m² (typical)
   - Formation enthalpies: -1 to -5 eV/atom (oxides)
   - Check against literature when possible

### Testing Workflow

1. **Create test script in `examples/`:**
   ```bash
   mkdir -p examples/my_feature
   vim examples/my_feature/test_my_feature.py
   ```

2. **Submit and wait:**
   ```bash
   python examples/my_feature/test_my_feature.py
   # Note the PK

   # Wait for completion (calculations take time)
   sleep 300  # 5 minutes, adjust as needed
   ```

3. **Check status:**
   ```bash
   verdi process show <PK>
   verdi process report <PK>
   ```

4. **If errors, fix and retry:**
   ```bash
   # Kill stuck processes
   verdi process kill <PK>

   # Kill VASP on current computer when testing (if needed)
   killall vasp_std

   # Restart daemon
   verdi daemon restart

   # Try again
   python examples/my_feature/test_my_feature.py
   ```

---

## 10. Troubleshooting

### Common Issues & Solutions

#### "Daemon not running"
```bash
verdi daemon start
verdi daemon status
```

#### "Code changes not taking effect"
```bash
verdi daemon restart
```

#### "Process stuck in Waiting"
```bash
# Check if sub-processes are running
verdi process report <PK>

# Check daemon logs
verdi daemon logshow
```

#### "VASP calculation failed"
```bash
# Check VASP output
verdi process show <VASP_CALC_PK>
verdi calcjob outputcat <VASP_CALC_PK>  # Shows stdout
verdi calcjob res <VASP_CALC_PK>        # Shows parsed results
```

#### "Unable to find total energy"
Check that VASP completed successfully:
```bash
verdi calcjob outputcat <PK>
# Look for "Total CPU time used" at the end
```

#### "Process excepted"
```bash
verdi process report <PK>
# Shows the exception traceback
```

#### "Cluster02 job not finishing"
```bash
# Kill VASP on remote
ssh cluster02 killall vasp_std

# Then kill the AiiDA process
verdi process kill <PK>
```

### Debugging Tips

1. **Start with `max_concurrent_jobs=1`** for debugging
2. **Check `verdi daemon logshow`** for daemon errors
3. **Use `verdi process report`** for detailed hierarchy
4. **Export structures** to verify they're correct:
   ```bash
   verdi data core.structure export <PK> -o check.cif
   ```

---

## 11. Code Style & Conventions

### Import Order

```python
# 1. Standard library
from __future__ import annotations
import typing as t
from copy import deepcopy

# 2. AiiDA
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida.engine import calcfunction

# 3. AiiDA-WorkGraph
from aiida_workgraph import task, WorkGraph, dynamic, namespace

# 4. Third-party scientific
import numpy as np
from pymatgen.core.surface import SlabGenerator
from ase.io import read

# 5. PS-TEROS
from teros.core.slabs import deep_merge_dicts
from teros.core.helper_functions import get_structure_from_file
```

### Type Hints

```python
# Use modern Python type hints
def my_function(
    structures: dict[str, orm.StructureData],
    parameters: orm.Dict,
    indices: list[int],
) -> orm.Float:

# For AiiDA-WorkGraph dynamic types
def scatter_function(
    inputs: t.Annotated[dict, dynamic(orm.StructureData)],
) -> t.Annotated[dict, namespace(
    outputs=dynamic(orm.Float),
)]:
```

### Docstrings

```python
def my_function(arg1: str, arg2: int = 10) -> orm.Dict:
    """
    Short description (one line).

    Longer description if needed, explaining the purpose
    and behavior of the function.

    Args:
        arg1: Description of arg1
        arg2: Description of arg2 (default: 10)

    Returns:
        Description of return value

    Raises:
        ValueError: When invalid input is provided

    Example:
        >>> result = my_function('test', arg2=20)
        >>> print(result.get_dict())
        {'key': 'value'}
    """
```

### Naming Conventions

| Type | Convention | Example |
|------|------------|---------|
| Modules | `snake_case` | `surface_energy.py` |
| Functions | `snake_case` | `calculate_surface_energy()` |
| Classes | `PascalCase` | `SlabGenerator` |
| Constants | `UPPER_CASE` | `WORKFLOW_PRESETS` |
| Private | `_prefix` | `_prepare_inputs()` |

### File Size Guidelines

- Keep modules under ~1000 lines
- Split large modules into submodules
- Each file should have a single responsibility

---

## Quick Reference Card

```bash
# Environment
source ~/envs/aiida/bin/activate
verdi profile set-default presto
verdi daemon restart

# Submit workflow
python examples/my_script.py

# Monitor
verdi process show <PK>
verdi process report <PK>

# Debug
verdi daemon logshow
verdi calcjob outputcat <PK>

# Control
verdi process kill <PK>
killall vasp_std
```
