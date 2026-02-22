# Comprehensive Guide: Creating Multi-Stage Run Scripts with the Lego Module

**Version:** 1.0  
**Date:** 2026-02-07  
**Target Audience:** AI agents and developers creating PS-TEROS workflows

---

## Table of Contents

1. [Overview](#1-overview)
2. [Architecture & Core Concepts](#2-architecture--core-concepts)
3. [Step-by-Step Guide](#3-step-by-step-guide)
4. [Brick Configuration Reference](#4-brick-configuration-reference)
5. [Common Patterns & Best Practices](#5-common-patterns--best-practices)
6. [Troubleshooting Guide](#6-troubleshooting-guide)
7. [Example Templates](#7-example-templates)

---

## 1. Overview

### 1.1 What is the Lego Module?

The Lego module (`teros.core.lego`) provides a **modular, brick-based system** for composing multi-stage computational workflows. Each "brick" represents a specific computational stage type (VASP relaxation, DOS calculation, Bader analysis, etc.), and bricks connect through a **typed port system** that ensures data flows correctly between stages.

### 1.2 Key Principles

1. **Explicit connections**: Every input declares where its data comes from (`structure_from`, `charge_from`, `restart`)
2. **Fail early with helpful errors**: Connection validation runs before AiiDA submission
3. **Pure Python metadata**: Port declarations in `connections.py` enable fast tier1 testing
4. **No magic**: Data flow is explicit, predictable, and debuggable

### 1.3 Available Brick Types

| Brick Type | Purpose | Produces Structure? | Key Use Cases |
|------------|---------|---------------------|---------------|
| `vasp` | Standard VASP calculations | ✅ (if NSW > 0) | Relaxation, SCF, static calculations |
| `dos` | Density of states | ❌ | Electronic structure analysis |
| `batch` | Parallel VASP with varying params | ❌ | Fukui functions, parameter scans |
| `bader` | Bader charge analysis | ❌ | Charge partitioning |
| `convergence` | ENCUT & k-points testing | ❌ | Finding optimal parameters |
| `thickness` | Slab thickness convergence | ❌ | Surface energy convergence |
| `hubbard_response` | Hubbard U response calculations | ❌ | DFT+U parameter determination |
| `hubbard_analysis` | Hubbard U regression | ❌ | Extract U from response data |
| `aimd` | Ab initio molecular dynamics | ✅ | MD simulations |
| `qe` | Quantum ESPRESSO calculations | ✅ | Alternative to VASP |
| `cp2k` | CP2K calculations | ✅ | Alternative to VASP |
| `generate_neb_images` | NEB image generation | ❌ | Create interpolated structures |
| `neb` | NEB pathway calculations | ✅ | Reaction pathways |

---

## 2. Architecture & Core Concepts

### 2.1 The Port System

Every brick declares **inputs** and **outputs** via a `PORTS` dictionary defined in `connections.py`:

```python
VASP_PORTS = {
    'inputs': {
        'structure': {
            'type': 'structure',
            'required': True,
            'source': 'auto',  # How to resolve this input
            'description': 'Atomic structure',
        },
        'restart_folder': {
            'type': 'remote_folder',
            'required': False,
            'source': 'restart',
            'description': 'Remote folder for WAVECAR/CHGCAR restart',
        },
    },
    'outputs': {
        'structure': {
            'type': 'structure',
            'conditional': {'incar_key': 'nsw', 'operator': '>', 'value': 0},
            'description': 'Relaxed structure (only if nsw > 0)',
        },
        'energy': {'type': 'energy', 'description': 'Total energy (eV)'},
        'misc': {'type': 'misc', 'description': 'Parsed VASP results'},
        'remote_folder': {'type': 'remote_folder'},
        'retrieved': {'type': 'retrieved'},
    },
}
```

**Port Types:**
- `structure`: Atomic structure (StructureData)
- `energy`: Total energy (Float)
- `misc`: Parsed calculation results (Dict)
- `remote_folder`: RemoteData link to cluster files
- `retrieved`: FolderData with retrieved files
- `dos_data`, `projectors`: DOS-specific outputs
- `bader_charges`: Bader analysis results
- `trajectory`: AIMD trajectory
- `convergence`: Convergence test results
- `hubbard_*`: Hubbard U workflow outputs
- `neb_images`: NEB image structures

### 2.2 Connection Resolution Modes

| Mode | Used By | Description |
|------|---------|-------------|
| `'auto'` | VASP | First stage uses input structure, later stages read `structure_from` |
| `'structure_from'` | DOS, batch, convergence | Reads `stage['structure_from']` |
| `'charge_from'` | Bader | Reads `stage['charge_from']` for CHGCAR/AECCAR files |
| `'restart'` | VASP, AIMD, NEB | Reads `stage['restart']` for RemoteData folder |
| `'ground_state_from'` | hubbard_response | Reads `stage['ground_state_from']` |
| `'response_from'` | hubbard_analysis | Reads `stage['response_from']` |

### 2.3 Stage Definition Structure

Every stage follows a **three-section pattern**:

```python
{
    # ═══ Section 1: Identity ═══
    'name': 'stage_name',       # Unique identifier
    'type': 'vasp',             # Brick type
    
    # ═══ Section 2: Connections ═══
    'structure_from': 'relax',  # Where to get structure (if applicable)
    'restart': None,            # Which stage to restart from (or None)
    
    # ═══ Section 3: Configuration ═══
    'incar': {...},             # Brick-specific parameters
    'kpoints_spacing': 0.03,
    'retrieve': ['CHGCAR'],
}
```

---

## 3. Step-by-Step Guide

### Step 1: Import Required Modules

```python
from pathlib import Path
from aiida import load_profile, orm
from pymatgen.core import Structure  # or ase.io import read
from teros.core.lego import quick_vasp_sequential, print_sequential_results
```

### Step 2: Load AiiDA Profile

```python
load_profile()  # Uses default profile, or load_profile('presto')
```

### Step 3: Load Structure(s)

**Option A: From file (Pymatgen)**
```python
struct_path = Path(__file__).parent / 'structure.vasp'
pmg_struct = Structure.from_file(str(struct_path))
structure = orm.StructureData(pymatgen=pmg_struct)
```

**Option B: From file (ASE)**
```python
from ase.io import read
structure = orm.StructureData(ase=read('structure.vasp'))
```

**Option C: Programmatically (Pymatgen)**
```python
from pymatgen.core import Lattice, Structure
lattice = Lattice.cubic(5.0)
struct = Structure(lattice, ['Fe'], [[0, 0, 0]])
structure = orm.StructureData(pymatgen=struct)
```

### Step 4: Define Cluster Configuration

**Local testing (localwork):**
```python
code_label = 'VASP-6.5.1@localwork'
potential_family = 'PBE'
potential_mapping = {'Sn': 'Sn_d', 'O': 'O'}
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 8,
    },
}
max_concurrent_jobs = 1  # CRITICAL: localwork runs ONE job at a time
```

**Production cluster (obelix):**
```python
code_label = 'VASP-6.5.1-idefix-4@obelix'
potential_family = 'PBE'
potential_mapping = {'Sn': 'Sn_d', 'O': 'O'}
options = {
    'resources': {
        'num_machines': 1,
        'num_mpiprocs_per_machine': 4,  # Hybrid MPI+OpenMP
    },
    'custom_scheduler_commands': '''#PBS -l cput=90000:00:00
#PBS -l nodes=1:ppn=88:skylake
#PBS -j oe
#PBS -N my_job_name''',
}
max_concurrent_jobs = 4  # Obelix can handle multiple parallel jobs
```

### Step 5: Define Stages

**Pattern: Define base INCAR, then stages**
```python
base_incar = {
    'encut': 520,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    'prec': 'Accurate',
    'algo': 'Normal',
}

stages = [
    # Stage 1: Relaxation
    {
        'name': 'relax',
        'type': 'vasp',
        'incar': {
            **base_incar,
            'ibrion': 2,
            'nsw': 100,
            'isif': 3,
            'lwave': False,
            'lcharg': False,
        },
        'restart': None,
        'kpoints_spacing': 0.03,
        'retrieve': ['CONTCAR', 'OUTCAR'],
    },
    # Stage 2: Static SCF on relaxed structure
    {
        'name': 'scf',
        'type': 'vasp',
        'structure_from': 'relax',
        'incar': {
            **base_incar,
            'nsw': 0,
            'ibrion': -1,
            'lwave': True,
            'lcharg': True,
        },
        'restart': None,
        'kpoints_spacing': 0.02,
        'retrieve': ['CHGCAR'],
    },
]
```

### Step 6: Submit Workflow

```python
result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label=code_label,
    potential_family=potential_family,
    potential_mapping=potential_mapping,
    options=options,
    max_concurrent_jobs=max_concurrent_jobs,
    name='my_workflow',
)

pk = result['__workgraph_pk__']
print(f"Submitted WorkGraph PK: {pk}")
print(f"Stages: {result['__stage_names__']}")
```

### Step 7: Monitor & Retrieve Results

**Monitor:**
```bash
verdi process show <PK>
verdi process report <PK>
verdi daemon logshow  # If something goes wrong
```

**Retrieve results (after completion):**
```python
from teros.core.lego import print_sequential_results, get_stage_results

# Print formatted summary
print_sequential_results(result)

# Extract specific stage results
scf_results = get_stage_results(result, 'scf')
print(f"SCF Energy: {scf_results['energy']:.6f} eV")
```

---

## 4. Brick Configuration Reference

### 4.1 VASP Brick

**Produces:** structure (if NSW > 0), energy, misc, remote_folder, retrieved

**Minimal configuration:**
```python
{
    'name': 'relax',
    'type': 'vasp',
    'incar': {'nsw': 100, 'ibrion': 2, 'isif': 3},
    'restart': None,
}
```

**Full configuration:**
```python
{
    'name': 'relax',
    'type': 'vasp',
    'structure': structure,      # Optional: override initial structure
    'structure_from': 'previous', # 'input', 'previous', or stage name
    'incar': {
        'encut': 520,
        'ediff': 1e-6,
        'ismear': 0,
        'sigma': 0.05,
        'ibrion': 2,
        'nsw': 100,
        'isif': 3,
        'prec': 'Accurate',
        'lreal': 'Auto',
        'lwave': False,
        'lcharg': False,
    },
    'restart': None,              # None or stage name for WAVECAR/CHGCAR restart
    'kpoints_spacing': 0.03,
    'fix_type': None,             # 'bottom', 'center', 'top' for selective dynamics
    'fix_thickness': 0.0,         # Thickness in Å (required if fix_type set)
    'fix_components': 'XYZ',      # 'XYZ', 'XY', 'Z'
    'supercell': None,            # [nx, ny, nz] for supercell transformation
    'retrieve': ['CONTCAR'],      # Additional files to retrieve
    'code_label': None,           # Override default code for this stage
    'options': None,              # Override scheduler options for this stage
}
```

**Key INCAR parameters:**
- `nsw`: 0 (static), 100+ (relaxation)
- `ibrion`: -1 (static), 2 (conjugate gradient), 1 (RMM-DIIS)
- `isif`: 2 (ions only), 3 (ions + cell)
- `ismear`: 0 (Gaussian, insulators), 1 (Methfessel-Paxton, metals)
- `lwave`: True (save WAVECAR), False (don't save)
- `lcharg`: True (save CHGCAR), False (don't save)

### 4.2 DOS Brick

**Produces:** energy (from SCF), scf_misc, dos_misc, dos_data, projectors

**Configuration:**
```python
{
    'name': 'dos',
    'type': 'dos',
    'structure_from': 'relax',    # REQUIRED
    'scf_incar': {
        'encut': 520,
        'ediff': 1e-6,
        'ismear': 0,
        'sigma': 0.05,
        'prec': 'Accurate',
        'algo': 'Normal',
    },
    'dos_incar': {
        'nedos': 2000,
        'lorbit': 11,              # Required for projected DOS
        'ismear': -5,              # Tetrahedron method recommended
        'prec': 'Accurate',
    },
    'kpoints_spacing': 0.03,       # SCF k-points
    'dos_kpoints_spacing': 0.02,   # DOS k-points (denser)
    'retrieve': ['DOSCAR'],
}
```

**Important:** BandsWorkChain automatically sets `lwave=True, lcharg=True` (SCF) and `icharg=11, istart=1` (DOS). Use lowercase INCAR keys.

### 4.3 Batch Brick

**Produces:** Multiple calculation results

**Configuration:**
```python
{
    'name': 'parameter_scan',
    'type': 'batch',
    'structure_from': 'relax',
    'calculations': {
        'encut_300': {
            'incar': {'encut': 300, 'nsw': 0},
            'kpoints_spacing': 0.05,
        },
        'encut_400': {
            'incar': {'encut': 400, 'nsw': 0},
        },
    },
    'base_incar': {'ediff': 1e-6, 'ismear': 0, 'sigma': 0.05},
    'retrieve': ['OUTCAR'],
}
```

### 4.4 Bader Brick

**Produces:** bader_charges

**Configuration:**
```python
{
    'name': 'bader',
    'type': 'bader',
    'charge_from': 'scf',
}
```

**Prerequisites:** The referenced stage MUST have `laechg: True` in INCAR and retrieve `['AECCAR0', 'AECCAR2', 'CHGCAR']`.

### 4.5 Hubbard Response & Analysis Bricks

**Response brick:**
```python
{
    'name': 'response',
    'type': 'hubbard_response',
    'ground_state_from': 'ground_state',
    'structure_from': 'relax',
    'target_species': 'Fe',
    'potential_values': [-0.2, -0.1, 0.1, 0.2],  # Do NOT include 0.0
    'ldaul': 2,  # 2=d-electrons, 3=f-electrons
    'ldauj': 0.0,
    'incar': {'encut': 520, 'ediff': 1e-6},
}
```

**Analysis brick:**
```python
{
    'name': 'analysis',
    'type': 'hubbard_analysis',
    'response_from': 'response',
    'structure_from': 'relax',
    'target_species': 'Fe',
    'ldaul': 2,
}
```

**Ground state requirements:**
```python
{
    'name': 'ground_state',
    'type': 'vasp',
    'incar': {
        'ldau': True,         # CRITICAL
        'ldautype': 3,
        'ldaul': [2, -1],
        'ldauj': [0.0, 0.0],
        'ldauu': [0.0, 0.0],  # U=0 for ground state
        'lmaxmix': 4,
        'lorbit': 11,
        'lwave': True,
        'lcharg': True,
    },
}
```

### 4.6 AIMD Brick

**Produces:** structure, energy, trajectory

**Configuration:**
```python
{
    'name': 'md',
    'type': 'aimd',
    'tebeg': 300,          # Target temperature (K)
    'nsw': 100,            # Number of MD steps
    'potim': 2.0,          # Time step (fs)
    'mdalgo': 2,           # Nose-Hoover thermostat
    'smass': 0.0,          # Nose mass (0.0 = auto)
    'supercell': [2, 2, 1],
    'incar': {'encut': 400, 'ediff': 1e-5},
    'restart': None,       # Or previous AIMD stage for velocity continuation
}
```

### 4.7 Convergence & Thickness Bricks

**Convergence:**
```python
{
    'name': 'conv_test',
    'type': 'convergence',
    'incar': {'prec': 'Accurate', 'ismear': 0, 'nsw': 0},
    'conv_settings': {
        'cutoff_start': 300,
        'cutoff_stop': 500,
        'cutoff_step': 100,
        'kspacing_start': 0.05,
        'kspacing_stop': 0.03,
        'kspacing_step': 0.01,
        'cutoff_kconv': 400,
        'kspacing_cutconv': 0.04,
    },
}
```

**Thickness (from previous VASP stage):**
```python
{
    'name': 'thick_conv',
    'type': 'thickness',
    'structure_from': 'bulk',
    'energy_from': 'bulk',
    'miller_indices': [1, 1, 0],
    'layer_counts': [3, 5, 7],
    'slab_incar': {'ibrion': 2, 'nsw': 50, 'isif': 2},
}
```

### 4.8 NEB Bricks

**Image generation:**
```python
{
    'name': 'make_images',
    'type': 'generate_neb_images',
    'initial_from': 'relax_initial',
    'final_from': 'relax_final',
    'n_images': 5,
    'method': 'idpp',
    'mic': True,
}
```

**NEB calculation:**
```python
{
    'name': 'neb_stage_1',
    'type': 'neb',
    'initial_from': 'relax_initial',
    'final_from': 'relax_final',
    'images_from': 'make_images',
    'incar': {'ibrion': 3, 'nsw': 150, 'iopt': 3, 'spring': -5, 'lclimb': False},
    'restart': None,
}
```

---

## 5. Common Patterns & Best Practices

### 5.1 Base INCAR Pattern

Define once, reuse with spread operator:
```python
base_incar = {
    'encut': 520,
    'ediff': 1e-6,
    'ismear': 0,
    'sigma': 0.05,
    'prec': 'Accurate',
}

stages = [
    {'name': 'relax', 'type': 'vasp', 'incar': {**base_incar, 'ibrion': 2, 'nsw': 100}},
    {'name': 'scf', 'type': 'vasp', 'structure_from': 'relax', 'incar': {**base_incar, 'nsw': 0}},
]
```

### 5.2 Testing Strategy

**Phase 1: Local testing (fast iteration)**
```python
code_label = 'VASP-6.5.1@localwork'
max_concurrent_jobs = 1
incar = {'encut': 300, 'ediff': 1e-4, 'nsw': 20}
kpoints_spacing = 0.06
```

**Phase 2: Production (validated parameters)**
```python
code_label = 'VASP-6.5.1-idefix-4@obelix'
max_concurrent_jobs = 4
incar = {'encut': 520, 'ediff': 1e-6, 'nsw': 100}
kpoints_spacing = 0.03
```

### 5.3 Common Workflows

**Relaxation → SCF → DOS:**
```python
stages = [
    {'name': 'relax', 'type': 'vasp', 'incar': {'ibrion': 2, 'nsw': 100}, 'restart': None},
    {'name': 'scf', 'type': 'vasp', 'structure_from': 'relax', 'incar': {'nsw': 0}, 'restart': None},
    {'name': 'dos', 'type': 'dos', 'structure_from': 'relax', 'scf_incar': {...}, 'dos_incar': {...}},
]
```

**Fukui function:**
```python
stages = [
    {'name': 'scf', 'type': 'vasp', 'incar': {'nsw': 0, 'lcharg': True}, 'restart': None},
    {
        'name': 'fukui',
        'type': 'batch',
        'structure_from': 'scf',
        'calculations': {
            'neutral': {'incar': {'nelect': 192}},
            'minus': {'incar': {'nelect': 193}},
        },
    },
]
```

**Bader analysis:**
```python
stages = [
    {
        'name': 'scf',
        'type': 'vasp',
        'incar': {'nsw': 0, 'lcharg': True, 'laechg': True},
        'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR'],
        'restart': None,
    },
    {'name': 'bader', 'type': 'bader', 'charge_from': 'scf'},
]
```

---

## 6. Troubleshooting Guide

### 6.1 Connection Errors

**Error: "structure_from references a 'dos' stage, which doesn't produce a structure"**

**Solution:** Only reference VASP, AIMD, QE, CP2K, or NEB stages for structures.

---

**Error: "charge_from='scf' requires laechg: True"**

**Solution:** Add required INCAR and retrieval:
```python
{'name': 'scf', 'type': 'vasp', 'incar': {'laechg': True}, 'retrieve': ['AECCAR0', 'AECCAR2', 'CHGCAR']}
```

---

**Error: "structure_from='relax' references a static calculation (nsw=0)"**

**Solution:** Ensure relaxation stage has `nsw > 0` or reference the input structure.

### 6.2 VASP Errors

**Error: "ICHARG=11 requires existing WAVECAR"**

**Solution:** Use DOS brick (handles automatically) or set `restart='scf'` with `lwave: True` in SCF stage.

---

**Error: "LDAU=True but LDAUU all zero"**

**Solution:** For Hubbard U ground state, set:
```python
'incar': {
    'ldau': True,
    'ldautype': 3,
    'ldaul': [2, -1],
    'ldauj': [0.0, 0.0],
    'ldauu': [0.0, 0.0],
    'lorbit': 11,
    'lwave': True,
    'lcharg': True,
}
```

### 6.3 Workflow Errors

**Error: "WorkGraph stuck in Waiting"**

**Solution:**
```bash
verdi process report <PK>  # Check dependencies
verdi daemon restart       # Restart daemon
```

**Error: "Code changes not taking effect"**

**Solution:**
```bash
verdi daemon restart  # ALWAYS after code changes
```

---

## 7. Example Templates

### 7.1 Simple Relaxation

```python
#!/usr/bin/env python
from pathlib import Path
from aiida import load_profile, orm
from pymatgen.core import Structure
from teros.core.lego import quick_vasp_sequential

load_profile()

structure = orm.StructureData(pymatgen=Structure.from_file('structure.vasp'))

stages = [
    {
        'name': 'relax',
        'type': 'vasp',
        'incar': {'encut': 520, 'ediff': 1e-6, 'ibrion': 2, 'nsw': 100, 'isif': 3},
        'restart': None,
        'kpoints_spacing': 0.03,
    },
]

result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label='VASP-6.5.1@localwork',
    potential_family='PBE',
    potential_mapping={'Fe': 'Fe_pv', 'O': 'O'},
    options={'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 8}},
    max_concurrent_jobs=1,
    name='simple_relax',
)

print(f"WorkGraph PK: {result['__workgraph_pk__']}")
```

### 7.2 Electronic Structure Pipeline

```python
#!/usr/bin/env python
from pathlib import Path
from aiida import load_profile, orm
from pymatgen.core import Structure
from teros.core.lego import quick_vasp_sequential

load_profile()

structure = orm.StructureData(pymatgen=Structure.from_file('structure.vasp'))

base_incar = {'encut': 520, 'ediff': 1e-6, 'ismear': 0, 'sigma': 0.05, 'prec': 'Accurate'}

stages = [
    {
        'name': 'relax',
        'type': 'vasp',
        'incar': {**base_incar, 'ibrion': 2, 'nsw': 100, 'isif': 3},
        'restart': None,
        'kpoints_spacing': 0.03,
    },
    {
        'name': 'scf',
        'type': 'vasp',
        'structure_from': 'relax',
        'incar': {**base_incar, 'nsw': 0, 'lwave': True, 'lcharg': True},
        'restart': None,
        'kpoints_spacing': 0.02,
    },
    {
        'name': 'dos',
        'type': 'dos',
        'structure_from': 'relax',
        'scf_incar': {**base_incar, 'nsw': 0},
        'dos_incar': {'encut': 520, 'prec': 'Accurate', 'nedos': 2000, 'lorbit': 11, 'ismear': -5},
        'kpoints_spacing': 0.03,
        'dos_kpoints_spacing': 0.02,
    },
]

result = quick_vasp_sequential(
    structure=structure,
    stages=stages,
    code_label='VASP-6.5.1-idefix-4@obelix',
    potential_family='PBE',
    potential_mapping={'Sn': 'Sn_d', 'O': 'O'},
    options={
        'resources': {'num_machines': 1, 'num_mpiprocs_per_machine': 4},
        'custom_scheduler_commands': '#PBS -l cput=90000:00:00\n#PBS -l nodes=1:ppn=88:skylake\n#PBS -j oe',
    },
    max_concurrent_jobs=4,
    name='electronic_structure',
)

print(f"WorkGraph PK: {result['__workgraph_pk__']}")
```

---

## 8. Validation Checklist

Before submitting, verify:

- [ ] All stage names are unique
- [ ] All `structure_from` references point to valid stages that produce structures
- [ ] All `restart` references point to valid stages
- [ ] Bader stages have upstream `laechg: True` and proper file retrieval
- [ ] Hubbard response has proper ground state (LDAU=True, LORBIT=11, LWAVE/LCHARG=True)
- [ ] AIMD stages have required parameters (tebeg, nsw, potim, mdalgo, smass)
- [ ] NEB stages have both endpoints and images source
- [ ] `max_concurrent_jobs` matches cluster (1 for localwork, 4+ for production)
- [ ] Code label matches cluster
- [ ] Potential family and mapping match available POTCARs

---

## 9. Additional Resources

### Documentation
- `teros/core/lego/README.md` - API reference
- `teros/core/lego/bricks/connections.py` - Port definitions
- `CLAUDE.md` / `AGENTS.md` - Development guides

### Examples
- `examples/lego/sno2_explorer/run_sequential_with_dos.py` - Relax + DOS
- `examples/lego/hubbard_u/run_sequential_relax_then_u.py` - Hubbard U
- `examples/lego/bader/run_bader_sno2.py` - Bader analysis
- `examples/lego/aimd/example_aimd.py` - AIMD
- `examples/lego/neb/run_neb_sno2.py` - NEB
- `examples/lego/thickness/run_thickness_sno2.py` - Thickness convergence

### Testing
```bash
pytest tests/test_lego_connections.py -m tier1 -v  # Connection tests
pytest tests/test_lego_*_integration.py -m tier2 -v  # Integration tests
```

---

**End of Guide**
