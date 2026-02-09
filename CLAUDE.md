# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Quick Reference

```bash
# Environment setup
source ~/envs/aiida/bin/activate
verdi profile set-default presto
verdi daemon restart  # CRITICAL: after any code changes

# Testing
pytest tests/ -v                              # Run all tests
pytest tests/test_imports.py -v               # Run specific test file
pytest tests/ -m tier1 -v                     # Fast unit tests (no AiiDA needed)
pytest tests/ -m "not slow" -v                # Exclude slow tests

# Linting (matches CI)
flake8 teros/ --max-line-length=120 --ignore=E501,W503,E402,F401

# Workflow monitoring
verdi process show <PK>                       # Show workflow status
verdi process report <PK>                     # Detailed hierarchy
verdi daemon logshow                          # Debug daemon issues
```

## Overview

**PS-TEROS** (Predicting Stability of TERminations of Oxide Surfaces) calculates surface Gibbs free energy using ab initio atomistic thermodynamics with AiiDA-WorkGraph.

**Technology Stack:** AiiDA (workflow/provenance), AiiDA-WorkGraph (task composition), AiiDA-VASP (DFT), Pymatgen (structures), ASE (I/O)

**Main entry point:** `teros.core.workgraph.build_core_workgraph()`

## Core Architecture

### Task Decorators (aiida_workgraph)

```python
from aiida import orm
from aiida_workgraph import task

@task.calcfunction  # Single computation, stored in AiiDA database
def extract_energy(misc: orm.Dict) -> orm.Float:
    return orm.Float(misc.get_dict()['total_energies']['energy_extrapolated'])

@task.graph(outputs=['energy', 'structure'])  # Composite workflow
def my_workflow(structure: orm.StructureData) -> dict:
    vasp_task = VaspTask(structure=structure, ...)
    return {'energy': extract_energy(misc=vasp_task.outputs.misc)}
```

### Scatter-Gather Pattern (Parallel Execution)

```python
from aiida_workgraph import dynamic, namespace
import typing as t

@task.graph(outputs=['energies'])
def relax_slabs_scatter(
    slab_structures: t.Annotated[dict, dynamic(orm.StructureData)],
) -> t.Annotated[dict, namespace(energies=dynamic(orm.Float))]:
    """Process {'term_0': StructureData, ...} → {'term_0': Float, ...}"""
    outputs = {'energies': {}}
    for label, slab in slab_structures.items():
        vasp_task = VaspTask(structure=slab, name=f'relax_{label}', ...)
        outputs['energies'][label] = extract_energy(misc=vasp_task.outputs.misc)
    return outputs
```

### Wrapping External WorkChains

```python
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task

VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
VaspTask = task(VaspWorkChain)  # Now usable in WorkGraphs
```

### Builder Input Pattern

All VASP calculations use a standardized `builder_inputs` dictionary:

```python
builder_inputs = {
    'parameters': {
        'incar': {
            'ENCUT': 520, 'EDIFF': 1e-6, 'ISMEAR': 0, 'SIGMA': 0.05,
            'IBRION': 2, 'NSW': 100, 'ISIF': 3,  # 3=full relax, 2=ions only
            'PREC': 'Accurate', 'LREAL': 'Auto', 'LWAVE': False, 'LCHARG': False,
        }
    },
    'options': {
        'resources': {'num_machines': 1, 'num_cores_per_machine': 24},
        'queue_name': 'normal', 'walltime_seconds': 86400,
    },
    'kpoints_spacing': 0.03,  # Å⁻¹
    'potential_family': 'PBE.54',
    'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
    'clean_workdir': False,
}
```

### Deep Merge for Parameter Overrides

```python
from teros.core.utils import deep_merge_dicts

base = {'parameters': {'incar': {'ENCUT': 400, 'ISMEAR': 0}}}
override = {'parameters': {'incar': {'ENCUT': 520}}}
result = deep_merge_dicts(base, override)
# {'parameters': {'incar': {'ENCUT': 520, 'ISMEAR': 0}}}  # ISMEAR preserved
```

## Workflow Presets

```python
from teros.core import build_core_workgraph

wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # One-liner activation
    structures_dir='/path/to/structures',
    ...
)
```

**Available presets:**
- `surface_thermodynamics` - Full oxide surface energy workflow
- `surface_thermodynamics_unrelaxed` - Quick screening, no slab relaxation
- `bulk_only` - Bulk structure optimization only
- `formation_enthalpy_only` - Formation enthalpy calculation
- `cleavage_only` - Cleavage energy for complementary surfaces
- `relaxation_energy_only` - Relaxation energy analysis
- `electronic_structure_bulk_only` - DOS/bands for bulk
- `electronic_structure_slabs_only` - DOS/bands for slabs
- `electronic_structure_bulk_and_slabs` - DOS/bands for both
- `aimd_only` - AIMD simulation
- `adsorption_energy` - Adsorption energy calculations
- `comprehensive` - Everything enabled

Override preset defaults with individual flags:
```python
workflow_preset='surface_thermodynamics',
compute_cleavage=True,   # Enable optional feature
max_concurrent_jobs=4,   # Control parallelism
```

## Module Structure

```
teros/core/
├── workgraph.py              # Main orchestration (~2000 lines), build_core_workgraph()
├── workflow_presets.py       # Preset definitions
├── thermodynamics.py         # Surface energy γ(Δμ_O, Δμ_M) calculations
├── hf.py                     # Formation enthalpy
├── cleavage.py               # Cleavage energy (complementary pairs)
├── slabs.py                  # Slab generation & relaxation
├── adsorption_energy.py      # E_ads = E_complete - E_substrate - E_molecule
├── fixed_atoms.py            # Selective dynamics (atom fixing)
├── utils.py                  # deep_merge_dicts, logging, helpers
├── constants.py              # Physical constants, unit conversions
├── aimd/                     # Ab initio molecular dynamics submodule
├── custom_calculation/       # Generic VASP calculations
├── lego/                     # Lightweight, incremental VASP calculations
├── surface_energy/           # Metal/intermetallic surface energy
├── surface_hydroxylation/    # OH/vacancy studies
├── builders/                 # Pre-configured parameter sets (Ag2O, Ag3PO4, etc.)
├── u_calculation/            # Hubbard U parameter (linear response method)
│   ├── __init__.py
│   ├── workgraph.py          # build_u_calculation_workgraph()
│   ├── tasks.py              # extract_d_electron_occupation, gather_responses, etc.
│   │                         # NOTE: Fixed 2026-02-02 - chi/chi_0 labels + sign convention
│   └── utils.py              # INCAR preparation, LDAU arrays, linear regression
│                             # NOTE: Fixed 2026-02-02 - potential sign negation
└── lego/                     # Multi-stage sequential workflow builder
    ├── __init__.py            # Public API: quick_vasp, quick_dos, quick_hubbard_u, etc.
    ├── workgraph.py           # quick_vasp_sequential(), _validate_stages(), _prepare_builder_inputs
    ├── tasks.py               # extract_energy, compute_dynamics calcfunctions
    ├── results.py             # Result extraction from completed WorkGraphs
    ├── utils.py               # get_status helper, lego-specific helpers
    └── bricks/                # Pluggable stage types (modular brick system)
        ├── __init__.py        # BRICK_REGISTRY, get_brick_module(), resolve_structure_from
        ├── connections.py     # PORTS dicts, PORT_TYPES, validate_connections() (pure Python)
        ├── vasp.py            # VASP brick (relaxation, SCF, etc.)
        ├── dos.py             # DOS brick (BandsWorkChain wrapper)
        ├── batch.py           # Batch brick (parallel VASP with varying params)
        ├── bader.py           # Bader brick (charge analysis)
        ├── cohp.py            # COHP brick (bonding analysis via LOBSTER)
        ├── convergence.py     # Convergence brick (ENCUT/k-points testing)
        ├── thickness.py       # Thickness brick (slab thickness convergence)
        ├── hubbard_response.py # Hubbard U response calculations (NSCF + SCF)
        ├── hubbard_analysis.py # Hubbard U regression and summary
        ├── aimd.py            # Ab initio molecular dynamics
        ├── qe.py              # Quantum ESPRESSO (PwBaseWorkChain)
        ├── cp2k.py            # CP2K (Cp2kBaseWorkChain)
        ├── generate_neb_images.py # NEB image interpolation/generation
        └── neb.py             # NEB workflow wrapper (vasp.neb)
```

**Submodule pattern:**
```
module_name/
├── __init__.py         # Public exports
├── workgraph.py        # Main WorkGraph builder
├── tasks.py            # @task.calcfunction helpers
└── utils.py            # Pure Python utilities (optional)
```

## Lego Module

The Lego module (`teros.core.lego`) provides a modular, brick-based system for composing VASP workflows from reusable stage types.

### Philosophy

The lego system treats computational stages as **physical Lego bricks**. Each brick has **studs** (outputs) and **sockets** (inputs) with specific shapes -- you can only snap together pieces that fit. VASP is the "joker brick": it produces every output type (structure, energy, remote folder, retrieved files), so any brick can connect to it. Specialized bricks (DOS, bader, batch, convergence) produce only their domain-specific outputs and cannot serve as structure sources.

Three core principles:

1. **Explicit connections.** Every input declares where its data comes from (`structure_from`, `charge_from`, `restart`). No implicit or magical data flows.
2. **Fail early with helpful errors.** `validate_connections()` runs before any AiiDA submission. Incompatible wiring (e.g., DOS to bader) gives a clear error naming the problem and suggesting which stages can provide what you need.
3. **Pure Python metadata.** All port declarations and validation logic live in `connections.py` with zero AiiDA imports, enabling tier1 testing.

### Brick Interface Contract

Each brick module in `bricks/` exports a `PORTS` dict plus 5 functions:

```python
PORTS                    # From connections.py (e.g., VASP_PORTS as PORTS)
validate_stage()         # Brick-specific config validation
create_stage_tasks()     # Build WorkGraph tasks
expose_stage_outputs()   # Wire up WorkGraph outputs
get_stage_results()      # Extract results from completed WorkGraph
print_stage_results()    # Format results for display
```

### Available Brick Types

| Type | Module | Description |
|------|--------|-------------|
| `vasp` | `bricks/vasp.py` | Standard VASP calculation (relax, static SCF) |
| `dos` | `bricks/dos.py` | Density of states via BandsWorkChain |
| `batch` | `bricks/batch.py` | Parallel VASP runs with varying parameters |
| `bader` | `bricks/bader.py` | Bader charge analysis |
| `cohp` | `bricks/cohp.py` | COHP bonding analysis via LOBSTER |
| `convergence` | `bricks/convergence.py` | ENCUT and k-points convergence testing |
| `thickness` | `bricks/thickness.py` | Slab thickness convergence testing |
| `hubbard_response` | `bricks/hubbard_response.py` | Hubbard U response calculations (NSCF + SCF per potential) |
| `hubbard_analysis` | `bricks/hubbard_analysis.py` | Hubbard U linear regression and summary |
| `aimd` | `bricks/aimd.py` | Ab initio molecular dynamics (IBRION=0) |
| `qe` | `bricks/qe.py` | Quantum ESPRESSO pw.x calculations (SCF, relax, vc-relax) |
| `cp2k` | `bricks/cp2k.py` | CP2K calculations (ENERGY, GEO_OPT, CELL_OPT, MD) |
| `generate_neb_images` | `bricks/generate_neb_images.py` | Generate intermediate NEB images from relaxed VASP endpoints |
| `neb` | `bricks/neb.py` | Run `vasp.neb` with images from generator stage or local folder |

### Port System (`connections.py`)

Each brick declares a `PORTS` dict with typed inputs and outputs. The `connections.py` module is **pure Python** (no AiiDA dependency) so it can be imported in tier1 tests.

```python
from teros.core.lego.bricks.connections import (
    validate_connections,    # Validate all inter-stage connections
    get_brick_info,          # Inspect a brick's PORTS
    ALL_PORTS,               # Registry of all PORTS dicts
    PORT_TYPES,              # Set of recognized port type strings
)
```

**Port types:** `structure`, `energy`, `misc`, `remote_folder`, `retrieved`, `dos_data`, `projectors`, `bader_charges`, `cohp_data`, `trajectory`, `convergence`, `file`, `hubbard_responses`, `hubbard_occupation`, `hubbard_result`, `neb_images`

**Source resolution modes:**
- `'auto'` -- VASP structure: first stage uses initial, then `'previous'`/`'input'`/explicit stage name
- `'structure_from'` -- reads `stage['structure_from']` (DOS, batch, convergence, hubbard)
- `'charge_from'` -- reads `stage['charge_from']` (bader)
- `'cohp_from'` -- reads `stage['cohp_from']` (cohp)
- `'restart'` -- reads `stage['restart']` (VASP restart folder)
- `'ground_state_from'` -- reads `stage['ground_state_from']` (hubbard_response)
- `'response_from'` -- reads `stage['response_from']` (hubbard_analysis)

### Connection Validation

`validate_connections(stages)` runs automatically before submission (called from `_validate_stages()`). It checks:

1. **Port type compatibility** -- referenced stage must produce the needed type
2. **Brick compatibility** -- `compatible_bricks` constraint (e.g., bader only from vasp)
3. **Prerequisites** -- upstream INCAR/retrieve requirements (e.g., bader needs `laechg: True`)
4. **Conditional outputs** -- warns when referencing a static stage's structure (nsw=0)
5. **Stage ordering** -- forward references are rejected

### Key Functions

| Function | Purpose |
|----------|---------|
| `quick_vasp(...)` | Single VASP calculation |
| `quick_vasp_batch(...)` | Multiple parallel VASP calcs on different structures |
| `quick_vasp_sequential(...)` | Multi-stage pipeline with restart chaining |
| `quick_dos(...)` | Single DOS calculation (SCF + DOS) |
| `quick_dos_batch(...)` | Multiple parallel DOS calcs |
| `quick_hubbard_u(...)` | Hubbard U calculation (ground state + response + analysis) |
| `quick_aimd(...)` | AIMD simulation (molecular dynamics) |
| `quick_qe(...)` | Single QE pw.x calculation |
| `quick_qe_sequential(...)` | Multi-stage QE pipeline with restart chaining |
| `get_results(pk)` | Extract results from single calc |
| `get_sequential_results(result)` | Extract all stage results |
| `get_stage_results(result, name)` | Extract one stage's results |
| `print_sequential_results(result)` | Print formatted results |

### Quick Functions

Convenience functions for common workflows:

```python
from teros.core.lego import quick_vasp, quick_dos, quick_hubbard_u, quick_qe, quick_vasp_sequential

# Single VASP calculation
wg = quick_vasp(structure, code_label, incar={...})

# DOS calculation
wg = quick_dos(structure, code_label, incar={...})

# Hubbard U calculation
wg = quick_hubbard_u(
    structure=structure,
    code_label='VASP-6.5.1@localwork',
    target_species='Fe',
    incar={'encut': 520, 'ediff': 1e-6, 'ismear': 0, 'sigma': 0.05},
    potential_values=[-0.2, -0.1, 0.1, 0.2],  # Default
    ldaul=2,  # 2=d-electrons, 3=f-electrons
)

# Single QE calculation
pk = quick_qe(
    structure=structure,
    code_label='pw@localhost',
    parameters={'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 50}},
    pseudo_family='SSSP/1.3/PBE/efficiency',
    kpoints_spacing=0.15,
    options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 4}},
)
```

### Sequential Workflows

Chain multiple stages with `quick_vasp_sequential`:

```python
from teros.core.lego import quick_vasp_sequential

stages = [
    {
        'name': 'relax',
        'type': 'vasp',
        'incar': {'encut': 520, 'ibrion': 2, 'nsw': 100, 'isif': 3},
        'restart': None,
    },
    {
        'name': 'ground_state',
        'type': 'vasp',
        'structure_from': 'relax',
        'incar': {
            'encut': 520, 'ediff': 1e-6, 'ismear': 0, 'sigma': 0.05,
            'ldau': False, 'lmaxmix': 4, 'lorbit': 11,
            'lwave': True, 'lcharg': True,
        },
        'restart': None,
        'retrieve': ['OUTCAR'],
    },
    {
        'name': 'response',
        'type': 'hubbard_response',
        'ground_state_from': 'ground_state',
        'structure_from': 'relax',
        'target_species': 'Sn',
        'potential_values': [-0.2, -0.1, 0.1, 0.2],
        'ldaul': 2,
        'incar': {'encut': 520, 'ediff': 1e-6, 'ismear': 0, 'sigma': 0.05},
    },
    {
        'name': 'analysis',
        'type': 'hubbard_analysis',
        'response_from': 'response',
        'structure_from': 'relax',
        'target_species': 'Sn',
        'ldaul': 2,
    },
]

wg = quick_vasp_sequential(structure, code_label, stages=stages, ...)
```

### Output Namespace Naming

Sequential workflow outputs use indexed prefixes for ordered display in `verdi process show`:

```
Outputs                       PK     Type
---------------------------   -----  -------------
s01_relax_2x2_rough
    vasp
        structure             34781  StructureData
        misc                  34780  Dict
        remote                34778  RemoteData
        retrieved             34779  FolderData
        energy                34785  Float
s02_relax_2x2_fine
    vasp
        ...
s03_dos_2x2
    scf
        misc                  35207  Dict
        ...
    dos
        misc                  35225  Dict
        ...
s04_fukui_minus_calcs_2x2
    batch
        delta_005
            misc              35143  Dict
            ...
        neutral
            ...
```

**Format:** `s{index:02d}_{stage_name}` (e.g., `s01_relax_rough`, `s02_dos`)

- The `s` prefix ensures valid Python identifiers (required by AiiDA link labels)
- Zero-padded index (`01`, `02`, ...) ensures alphabetical sort matches stage execution order
- Stage name provides descriptive context

**Brick-specific output structure:**

| Brick Type | Output Namespace Structure |
|------------|---------------------------|
| `vasp` | `s{N}_{name}.vasp.{energy,structure,misc,remote,retrieved}` |
| `dos` | `s{N}_{name}.scf.{misc,remote,retrieved}` + `s{N}_{name}.dos.{misc,remote,retrieved}` |
| `batch` | `s{N}_{name}.batch.{calc_label}.{energy,misc,remote,retrieved}` |

### AIMD Brick Status (2026-02-06)

The AIMD VASP brick is working great now for multi-stage MD workflows.

- Trajectory positions are normalized to Cartesian coordinates before exposure/concatenation.
- Combined trajectory output uses indexed naming: `s{N+1:02d}_combined_trajectory`
  (example: after 6 AIMD stages, output is `s07_combined_trajectory`).
- Velocity continuation across AIMD restart stages is supported via POSCAR velocity injection.
- Viewer helper script available at `~/.local/bin/trajview`:
  `trajview <TrajectoryData_PK>` opens the trajectory in ASE view.

### Hubbard U Stage Configuration

The Hubbard U workflow uses two brick types that wrap `teros.core.u_calculation`:

**`hubbard_response` brick** - runs NSCF + SCF response calculations:

| Key | Required | Default | Description |
|-----|----------|---------|-------------|
| `ground_state_from` | Yes | — | Name of the ground state vasp stage |
| `structure_from` | Yes | — | `'input'` or name of previous stage |
| `target_species` | Yes | — | Element symbol (e.g., 'Fe', 'Sn') |
| `potential_values` | No | `[-0.2, -0.1, 0.1, 0.2]` | Perturbation potentials (eV), no 0.0 |
| `ldaul` | No | `2` | Angular momentum (2=d, 3=f) |
| `ldauj` | No | `0.0` | Exchange J parameter |
| `incar` | No | `{}` | Base INCAR parameters |

**`hubbard_analysis` brick** - linear regression and summary (no VASP):

| Key | Required | Default | Description |
|-----|----------|---------|-------------|
| `response_from` | Yes | — | Name of the hubbard_response stage |
| `structure_from` | Yes | — | `'input'` or name of previous stage |
| `target_species` | Yes | — | Element symbol (e.g., 'Fe', 'Sn') |
| `ldaul` | No | `2` | Angular momentum (2=d, 3=f) |

The ground state is a standard vasp brick with `lorbit: 11, lwave: True, lcharg: True, ldau: False`. The `quick_hubbard_u()` convenience function builds all 3 stages automatically.

### CP2K Stage Configuration

The CP2K brick wraps `Cp2kBaseWorkChain` from aiida-cp2k. Key differences from VASP:
- Parameters are nested dicts (GLOBAL, FORCE_EVAL, MOTION)
- Energy is in Hartree, auto-converted to eV
- Restart uses `parent_calc_folder` (RemoteData)
- Fixed atoms use `MOTION.CONSTRAINT.FIXED_ATOMS` with 1-based LIST

**`cp2k` brick** - CP2K calculations (ENERGY, GEO_OPT, CELL_OPT, MD):

| Key | Required | Default | Description |
|-----|----------|---------|-------------|
| `parameters` | Yes | — | Nested CP2K dict (GLOBAL, FORCE_EVAL, MOTION) |
| `restart` | Yes | — | None or previous CP2K stage name |
| `file` | Yes* | — | Dict with 'basis' and 'pseudo' keys |
| `basis_file` | Yes* | — | Path or SinglefileData for BASIS_MOLOPT |
| `pseudo_file` | Yes* | — | Path or SinglefileData for GTH_POTENTIALS |
| `structure_from` | No | 'previous' | 'input', 'previous', or stage name |
| `fix_type` | No | None | 'bottom', 'center', 'top' |
| `fix_thickness` | No | 0.0 | Thickness in Å (required if fix_type set) |
| `fix_components` | No | 'XYZ' | Components to fix ('XYZ', 'XY', 'Z') |
| `supercell` | No | None | [nx, ny, nz] for supercell transformation |
| `max_iterations` | No | 3 | Max workchain iterations |
| `retrieve` | No | [] | Additional files to retrieve |
| `code_label` | No | context | Override default code for this stage |
| `options` | No | context | Override scheduler options for this stage |

*Either `file` dict or both `basis_file` + `pseudo_file` required.

**Conditional outputs:**
- `structure` only available when `RUN_TYPE` is GEO_OPT, CELL_OPT, or MD
- `trajectory` only available when `RUN_TYPE` is MD

## VASP Integration

### VaspWorkChain Outputs

| Output | Type | Description |
|--------|------|-------------|
| `structure` | `StructureData` | Relaxed structure (if NSW > 0) |
| `misc` | `Dict` | Parsed results (energies, forces) |
| `remote_folder` | `RemoteData` | Link to files on cluster |

### Accessing Energy

```python
misc_dict = vasp_task.outputs.misc.get_dict()
energy = misc_dict['total_energies']['energy_extrapolated']  # Recommended
```

### Key INCAR Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `ISMEAR` | Smearing | 0 (insulators), 1 (metals) |
| `ISIF` | Stress control | 2 (ions), 3 (ions+cell) |
| `IBRION` | Relaxation | 2 (CG), -1 (static) |
| `NSW` | Max ionic steps | 0 (SCF), 100+ (relax) |

### K-points Spacing Convention

| Accuracy Level | `kpoints_spacing` (Å⁻¹) |
|----------------|-------------------------|
| Gamma-Only     | 10                      |
| Low            | 0.06 - 0.04             |
| Medium         | 0.04 - 0.03             |
| Fine           | 0.02 - 0.01             |

**0.03 - 0.04 is generally precise enough for most calculations.**

## Testing

### Testing Environment by Location

**IMPORTANT:** Choose the testing environment based on the current working directory:

| Working Directory | Testing Strategy | Code Label |
|-------------------|------------------|------------|
| `/home/thiagotd` (home computer) | Jump directly to production on obelix | `VASP-6.5.1-idefix-4@obelix` |
| `/home/trevizam` (work computer) | 1. Test locally first, then production | `vasp-6.5.1-std@localhost` → `VASP-6.5.1-idefix-4@obelix` |

**Obelix cluster configuration:**
```python
code_label = 'VASP-6.5.1-idefix-4@obelix'
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 4,  # PROCESS_MPI=4 (hybrid MPI+OpenMP)
    },
    'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N MyJobName''',
}
potential_family = 'PBE'
```

### Test Markers

```python
@pytest.mark.tier1        # Fast unit tests (no AiiDA)
@pytest.mark.tier2        # Integration tests (AiiDA mock)
@pytest.mark.tier3        # End-to-end (real calculations)
@pytest.mark.slow         # Long-running tests
@pytest.mark.requires_aiida
```

### Test Cluster (localwork)

```python
code_label = 'VASP-6.5.1@localwork'
max_concurrent_jobs = 1  # CRITICAL: localwork runs ONE job at a time
```

### Validation Checklist

1. Process returns `Finished [0]}` via `verdi process show <PK>`
2. Expected outputs present in "Outputs" section
3. Values physically reasonable (surface energies: 0.5-3.0 J/m²)

## AiiDA Commands

```bash
# Monitoring
verdi process show <PK>           # Status
verdi process report <PK>         # Hierarchy
verdi process list -a -p 1        # Recent processes

# Control
verdi process kill <PK>           # Kill process
verdi daemon restart              # After code changes!

# Inspect
verdi calcjob outputcat <PK>      # VASP stdout
verdi data core.structure export <PK> -o out.cif
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Code changes not taking effect | `verdi daemon restart` |
| Process stuck in Waiting | `verdi process report <PK>` to check sub-processes |
| VASP calculation failed | `verdi calcjob outputcat <PK>` for stdout |
| Process excepted | `verdi process report <PK>` for traceback |

## Adding New Features

### Adding a Parameter to build_core_workgraph

1. Add to function signature with default: `my_param: bool = False`
2. Use conditionally: `if my_param: ...`
3. Document in docstring
4. Update CHANGE.md

### Creating a New Submodule

1. Create `teros/core/my_module/` with `__init__.py`, `workgraph.py`, `tasks.py`
2. Implement `@task.calcfunction` in tasks.py
3. Implement `@task.graph` and `build_*_workgraph()` in workgraph.py
4. Export in `__init__.py` and `teros/core/__init__.py`
5. Add example in `examples/`

### Adding a New Lego Brick

When adding a new brick type to the lego module, follow this procedure. Every new brick must ship with PORTS, tests, an example script, and updated documentation.

#### Step 1: Declare ports in `connections.py`

Add a `MYBRICK_PORTS` dict to `teros/core/lego/bricks/connections.py`. This file has **no AiiDA imports** — keep it that way.

```python
MYBRICK_PORTS = {
    'inputs': {
        'structure': {
            'type': 'structure',         # Must be in PORT_TYPES
            'required': True,
            'source': 'structure_from',  # Which stage config field to read
            'description': 'Input structure',
        },
        # Add more inputs as needed. Use 'compatible_bricks' and
        # 'prerequisites' if this brick requires specific upstream settings.
    },
    'outputs': {
        'my_result': {
            'type': 'misc',              # Or add a new type to PORT_TYPES
            'description': 'Main output',
        },
        # Only declare outputs the brick actually produces.
        # Do NOT declare 'structure' unless the brick genuinely
        # creates a new/modified structure (decision: no pass-through).
    },
}
```

Then register it in `ALL_PORTS`:

```python
ALL_PORTS = {
    'vasp': VASP_PORTS,
    'dos': DOS_PORTS,
    ...
    'mybrick': MYBRICK_PORTS,   # Add here
}
```

If the brick introduces a new output type (e.g., `'phonon_bands'`), add it to `PORT_TYPES` at the top of the file.

#### Step 2: Write unit tests for the ports and connections

Add tests to `tests/test_lego_connections.py` **before** implementing the brick logic. Tests import from `connections.py` directly (no AiiDA needed), so they run fast. Cover at minimum:

- Port type recognition: all input/output types are in `PORT_TYPES`
- Output count and names match expectations
- `validate_connections()` passes for a valid pipeline using the new brick
- `validate_connections()` rejects invalid connections (wrong source brick, missing required field, forward reference)
- If the brick has `prerequisites`, test that missing prerequisites raise errors
- If the brick has `compatible_bricks`, test that incompatible sources are rejected
- Run with: `pytest tests/test_lego_connections.py -v`

#### Step 3: Create the brick module

Create `teros/core/lego/bricks/mybrick.py` with the 5 required functions + PORTS import:

```python
from .connections import MYBRICK_PORTS as PORTS


def validate_stage(stage: dict, stage_names: set) -> None:
    """Validate brick-specific config (called per-stage before connections)."""
    ...

def create_stage_tasks(wg, stage, stage_name, context) -> dict:
    """Create WorkGraph tasks. Use resolve_structure_from() for structures."""
    from . import resolve_structure_from
    ...

def expose_stage_outputs(wg, stage_name, stage_tasks_result) -> None:
    """Wire outputs onto the WorkGraph using: wg.outputs.name = socket"""
    ...

def get_stage_results(wg_node, wg_pk, stage_name) -> dict:
    """Extract results from a completed WorkGraph node."""
    ...

def print_stage_results(index, stage_name, stage_result) -> None:
    """Print formatted results to stdout."""
    ...
```

#### Step 4: Register the brick

In `teros/core/lego/bricks/__init__.py`:

```python
from . import vasp, dos, batch, bader, convergence, mybrick

BRICK_REGISTRY = {
    ...
    'mybrick': mybrick,
}
```

Also update the import in `__init__.py` to re-export any new symbols from `connections.py` if needed.

#### Step 5: Create an example script

Add `examples/lego/example_mybrick.py` showing a real pipeline that uses the new brick. The example should be runnable (with the right VASP code configured) and show the typical stage configuration pattern.

#### Step 6: Run all tests and lint

```bash
pytest tests/test_lego_connections.py -v   # Port/connection tests (tier1, fast)
pytest tests/test_lego_bricks.py -v        # Brick validate_stage tests (tier1)
flake8 teros/core/lego/bricks/ --max-line-length=120 --ignore=E501,W503,E402,F401
```

#### Step 7: Update CLAUDE.md

Update the following sections in this file:

1. **Module Structure tree** — add the new `.py` file under `bricks/`
2. **Brick Types table** — add a row with the brick name, purpose, and key config fields
3. **Port types list** — if you added a new port type to `PORT_TYPES`

#### Checklist

- [ ] `MYBRICK_PORTS` in `connections.py` with correct types and sources
- [ ] Registered in `ALL_PORTS` dict
- [ ] New port types (if any) added to `PORT_TYPES`
- [ ] Unit tests in `test_lego_connections.py` (ports, valid pipelines, rejections)
- [ ] `mybrick.py` with 5 functions + PORTS import
- [ ] Registered in `BRICK_REGISTRY` in `__init__.py`
- [ ] Example script in `examples/lego/`
- [ ] All tests pass, lint clean
- [ ] CLAUDE.md updated (tree, brick table, port types)

### Module Development & Testing Procedure

Follow this standardized workflow when developing new modules:

#### Phase 1: Local Testing

1. **Create a minimal test script** in `examples/` using:
   - Simple test structure (e.g., primitive cell, 2-4 atoms)
   - Local VASP code: `VASP-6.5.1@localwork` (8 processors)
   - Light parameters: low ENCUT, coarse k-points, few ionic steps
   ```python
   # Example test configuration
   code_label = 'VASP-6.5.1@localwork'
   options = {
       'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8},
       # Note: use num_mpiprocs_per_machine, NOT num_cores_per_machine
   }
   potential_family = 'PBE'  # Check available families with verdi data core.potcar listfamilies
   ```

2. **Submit and monitor:**
   ```bash
   verdi daemon restart  # CRITICAL: always restart after code changes
   python examples/my_module/test_script.py
   verdi process show <PK>
   verdi process report <PK>  # For detailed task hierarchy
   ```

3. **Debug iteratively** until `Finished [0]}`:
   - Check `verdi calcjob outputcat <PK>` for VASP errors
   - Check `verdi process report <PK>` for Python tracebacks
   - Common issues (see troubleshooting table below)

#### Phase 2: Validate Results

1. **Compare with reference values** (literature, VASP wiki, known benchmarks)
2. **Check physical reasonableness:**
   - Surface energies: 0.5-3.0 J/m²
   - Formation enthalpies: negative for stable compounds
   - Band gaps: compare with experimental values

#### Phase 3: Expose Outputs

1. **Add WorkGraph outputs** using the correct pattern:
   ```python
   # CORRECT: Use direct assignment
   wg.outputs.summary = summary_task.outputs.result
   wg.outputs.energy = energy_task.outputs.result
   wg.outputs.structure = vasp_task.outputs.structure

   # WRONG: Don't use add_output (outputs won't appear in verdi process show)
   # wg.add_output('any', 'summary', from_socket=summary_task.outputs.result)
   ```

2. **Note on VaspWorkChain outputs:**
   - `structure` only available if NSW > 0 (relaxation)
   - Static SCF (NSW=0) only produces: `misc`, `remote_folder`, `retrieved`

3. **Create helper functions** for result extraction:
   ```python
   def get_my_results(workgraph) -> dict:
       """Extract results from completed WorkGraph (accepts PK, node, or live WorkGraph)."""
       from aiida.common.links import LinkType

       if isinstance(workgraph, (int, str)):
           workgraph = orm.load_node(workgraph)

       # For stored nodes, traverse CALL_CALC links to find calcfunction outputs
       if hasattr(workgraph, 'base'):
           calcfuncs = workgraph.base.links.get_outgoing(link_type=LinkType.CALL_CALC)
           for link in calcfuncs.all():
               if link.link_label == 'my_summary_task':
                   for out_link in link.node.base.links.get_outgoing(
                       link_type=LinkType.CREATE
                   ).all():
                       if out_link.link_label == 'result':
                           return out_link.node.get_dict()

       raise ValueError("Could not find results")

   def print_my_summary(workgraph) -> None:
       """Print formatted summary."""
       results = get_my_results(workgraph)
       # ... formatted output
   ```

#### Phase 4: Final Verification

1. **Run a fresh test** with all changes:
   ```bash
   verdi daemon restart
   python examples/my_module/test_script.py
   ```

2. **Verify outputs appear in `verdi process show <PK>`:**
   ```
   Outputs              PK  Type
   ----------------  -----  ------
   summary           12345  Dict
   energy            12346  Float
   structure         12347  StructureData
   ```

3. **Test helper functions:**
   ```python
   from teros.core.my_module import print_my_summary, get_my_results
   print_my_summary(<PK>)
   results = get_my_results(<PK>)
   ```

#### Optional: Production Testing on Cluster

> **Note:** Only perform production testing when explicitly requested. Local testing is sufficient for most development iterations.

If production testing is needed, use cluster code with PBS scheduler. **The resource key differs by cluster:**

**Obelix cluster** — uses `num_mpiprocs_per_machine`:
```python
code_label = 'VASP-6.5.1-idefix-4@obelix'
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 4,
    },
    'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N MyJobName''',
}
```

**Bohr cluster** — uses `num_cores_per_machine` (not `num_mpiprocs_per_machine`):
- **par40 queue**: 1 node, 40 processors
- **par120 queue**: 3 nodes, 40 processors per node
- **teste queue**: 1 node, 5 processors (testing only)
- **gpu_a100 queue**: 1 node, 56 processors + 1× A100 GPU

```python
code_label = 'VASP-VTST-6.4.3@bohr'
# par40 queue (1 node)
options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 40,
    },
    'custom_scheduler_commands': '#PBS -q par40\n#PBS -j oe\n#PBS -N MyJobName',
}
# par120 queue (3 nodes)
options = {
    'resources': {
        'num_machines': 3,
        'num_cores_per_machine': 40,
    },
    'custom_scheduler_commands': '#PBS -q par120\n#PBS -j oe\n#PBS -N MyJobName',
}
```

**Lovelace cluster** — uses `num_cores_per_machine` (same as bohr):
- **par128 queue**: 1 node, 128 processors
- **parexp queue**: 1–11 nodes, 48 processors per node
- **paralela queue**: 2–8 nodes, 128 processors per node (min 256 processors)
- **testes queue**: 1 node, 2–4 processors (testing only)
- **umagpu queue**: 1 node, 16 processors + 1 GPU
- **duasgpus queue**: 1 node, 32 processors + 2 GPUs

```python
code_label = 'VASP-6.5.0@lovelace'
# par128 queue (1 node)
options = {
    'resources': {
        'num_machines': 1,
        'num_cores_per_machine': 128,
    },
    'custom_scheduler_commands': '#PBS -q par128\n#PBS -j oe\n#PBS -N MyJobName',
}
# parexp queue (e.g. 4 nodes)
options = {
    'resources': {
        'num_machines': 4,
        'num_cores_per_machine': 48,
    },
    'custom_scheduler_commands': '#PBS -q parexp\n#PBS -j oe\n#PBS -N MyJobName',
}
# paralela queue (e.g. 2 nodes, min 2)
options = {
    'resources': {
        'num_machines': 2,
        'num_cores_per_machine': 128,
    },
    'custom_scheduler_commands': '#PBS -q paralela\n#PBS -j oe\n#PBS -N MyJobName',
}
```

**GAM codes:** AiiDA code labels containing `GAM` (e.g., `VASPGAM-6.5.0@lovelace`, `VASP-6.5.1-idefix-4-GAM@obelix`) use the `vasp_gam` binary instead of `vasp_std`. Use GAM codes for Gamma-only calculations (k-points = [1, 1, 1]) for significant speedup. Use `vasp_std` codes for any calculation with multiple k-points.

Use production-quality parameters (proper ENCUT, converged k-points, appropriate supercell) and verify results match expectations.

#### Common Issues During Development

| Issue | Cause | Solution |
|-------|-------|----------|
| INCAR keys not recognized | AiiDA-VASP requires lowercase | Use `'encut'` not `'ENCUT'` |
| POTCAR family not found | Wrong family name | Check with `verdi data core.potcar listfamilies` |
| Job requests too many CPUs | Wrong resource key | Cluster-dependent: bohr uses `num_cores_per_machine`; obelix/localwork use `num_mpiprocs_per_machine` |
| OUTCAR parsing fails | Parser doesn't extract needed data | Parse OUTCAR directly via `retrieved` FolderData |
| Outputs not in `verdi process show` | Used `wg.add_output()` | Use `wg.outputs.name = socket` instead |
| Results extraction fails on stored node | Accessing `.tasks` on AiiDA node | Traverse `CALL_CALC` links instead |
| V=0 baseline shift | Different INCAR settings at V=0 | Ensure consistent settings across all calculations |

## Code Style

**Import order:** stdlib → aiida → aiida_workgraph → third-party → teros

**Type hints:** Use `t.Annotated[dict, dynamic(orm.StructureData)]` for WorkGraph dynamic types

**Docstrings:** Google style with Args, Returns, Example sections

## CI/CD

CI runs on push to `main`/`develop` and PRs. Tests Python 3.9, 3.10, 3.11 with PostgreSQL and RabbitMQ.

```bash
# Local CI simulation
pip install -e ".[dev]"
python -m py_compile teros/core/*.py
pytest tests/ -v
flake8 teros/ --max-line-length=120 --ignore=E501,W503,E402,F401
```
