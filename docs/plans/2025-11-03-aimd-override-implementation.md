# AIMD Override System Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Enable per-structure INCAR parameter overrides in standalone AIMD module

**Architecture:** Modify `aimd_single_stage_scatter()` to accept base parameters plus optional per-structure overrides. Standalone module extracts INCAR overrides from three levels (structure/stage/matrix) and passes merged dict to scatter function.

**Tech Stack:** AiiDA, AiiDA-WorkGraph, Python 3.13, pytest

**Design:** `docs/plans/2025-11-03-aimd-override-system.md`

---

## Task 1: Update aimd_single_stage_scatter Signature

**Files:**
- Modify: `teros/core/aimd_functions.py:122-135`

**Step 1: Update function signature**

Find the `aimd_single_stage_scatter` function signature (starts around line 122) and modify:

```python
@task.graph
def aimd_single_stage_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    temperature: float,
    steps: int,
    code: orm.Code,
    base_aimd_parameters: dict,  # RENAMED from aimd_parameters
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
    restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {},
    structure_aimd_overrides: dict[str, dict] = None,  # NEW PARAMETER
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(structures=dynamic(orm.StructureData), remote_folders=dynamic(orm.RemoteData), energies=dynamic(orm.Float))]:
```

**Step 2: Update docstring**

Update the docstring (lines ~136-161) to document new parameters:

```python
    """
    Run single AIMD stage on all slabs in parallel using scatter-gather pattern.

    This function handles ONE temperature/timestep stage for all slabs.
    Call it multiple times sequentially to build multi-stage AIMD workflows.

    Args:
        slabs: Dictionary of slab structures to run AIMD on
        temperature: Target temperature in K
        steps: Number of MD steps for this stage
        code: VASP code
        base_aimd_parameters: Base AIMD INCAR parameters (IBRION=0, MDALGO, etc.)
                             Applied to all structures by default
        potential_family: Potential family name
        potential_mapping: Element to potential mapping
        options: Scheduler options
        kpoints_spacing: K-points spacing
        clean_workdir: Whether to clean work directory
        restart_folders: Optional dict of RemoteData for restart (from previous stage)
        structure_aimd_overrides: Optional per-structure INCAR overrides.
                                 Format: {structure_name: {INCAR_key: value}}
                                 Missing structures use base_aimd_parameters.
        max_number_jobs: Maximum number of concurrent VASP calculations (None = unlimited)

    Returns:
        Dictionary with outputs per slab:
            - structures: Output structures from this stage
            - remote_folders: RemoteData nodes for potential next stage restart
            - energies: Total energies from this stage
    """
```

**Step 3: Commit signature change**

```bash
git add teros/core/aimd_functions.py
git commit -m "refactor(aimd): update aimd_single_stage_scatter signature for overrides

- Rename aimd_parameters -> base_aimd_parameters
- Add structure_aimd_overrides parameter
- Update docstring"
```

---

## Task 2: Implement Override Merging Logic

**Files:**
- Modify: `teros/core/aimd_functions.py:178-181`

**Step 1: Add merge logic in scatter loop**

Find the scatter loop (starts around line 178) and modify the parameter preparation:

```python
    # Scatter: create AIMD task for each slab (runs in parallel)
    for slab_label, slab_structure in slabs.items():
        # Merge base parameters with structure-specific overrides
        if structure_aimd_overrides and slab_label in structure_aimd_overrides:
            # Shallow merge: structure overrides take precedence
            merged_params = {**base_aimd_parameters, **structure_aimd_overrides[slab_label]}
        else:
            # No override for this structure, use base
            merged_params = base_aimd_parameters

        # Prepare parameters for this stage
        stage_params = prepare_aimd_parameters(merged_params, temperature, steps)
```

Replace the old line:
```python
        stage_params = prepare_aimd_parameters(aimd_parameters, temperature, steps)
```

**Step 2: Commit merge logic**

```bash
git add teros/core/aimd_functions.py
git commit -m "feat(aimd): implement per-structure INCAR override merging

Merge logic:
- If structure has override, shallow merge with base
- Missing structures use base parameters only
- Override values replace base for same INCAR key"
```

---

## Task 3: Update Main Workflow Caller

**Files:**
- Modify: `teros/core/workgraph.py:1660`

**Step 1: Find and update parameter name**

Find the import and call to `aimd_single_stage_scatter` (around line 1660):

```python
        # Import appropriate scatter function based on calculator
        if calculator == 'vasp':
            from teros.core.aimd_functions import aimd_single_stage_scatter
            aimd_scatter_func = aimd_single_stage_scatter
```

Find where it's called (search for `aimd_scatter_func(` or check AIMD integration section). Update the call to use new parameter name:

```python
        # Change this line:
        aimd_parameters=aimd_base_params,

        # To this:
        base_aimd_parameters=aimd_base_params,
        structure_aimd_overrides=None,
```

**Step 2: Commit main workflow update**

```bash
git add teros/core/workgraph.py
git commit -m "fix(workgraph): update aimd_single_stage_scatter call for new signature

- Rename aimd_parameters -> base_aimd_parameters
- Add structure_aimd_overrides=None parameter"
```

---

## Task 4: Implement Override Extraction in Standalone Module

**Files:**
- Modify: `teros/core/aimd/workgraph.py:175-230`

**Step 1: Replace stage loop with override extraction**

Find the stage loop in `build_aimd_workgraph` (starts around line 175). Replace the section from "for stage_idx, stage_config" through the `stage_task = wg.add_task(...)` call:

```python
    for stage_idx, stage_config in enumerate(aimd_stages):
        temperature = stage_config['temperature']
        steps = stage_config['steps']

        # Build per-structure INCAR overrides for this stage
        structure_incar_overrides = {}

        for struct_name in prepared_structures:
            # Start with empty override
            override = {}

            # Apply structure-level override (if exists)
            if structure_overrides and struct_name in structure_overrides:
                struct_override = structure_overrides[struct_name]
                if 'parameters' in struct_override and 'incar' in struct_override['parameters']:
                    override.update(struct_override['parameters']['incar'])

            # Apply stage-level override (if exists)
            if stage_overrides and stage_idx in stage_overrides:
                stage_override = stage_overrides[stage_idx]
                if 'parameters' in stage_override and 'incar' in stage_override['parameters']:
                    override.update(stage_override['parameters']['incar'])

            # Apply matrix-level override (highest priority)
            matrix_key = (struct_name, stage_idx)
            if matrix_overrides and matrix_key in matrix_overrides:
                matrix_override = matrix_overrides[matrix_key]
                if 'parameters' in matrix_override and 'incar' in matrix_override['parameters']:
                    override.update(matrix_override['parameters']['incar'])

            # Only add to dict if we have actual overrides
            if override:
                structure_incar_overrides[struct_name] = override

        # Extract base INCAR from builder_inputs
        if 'parameters' in builder_inputs and 'incar' in builder_inputs['parameters']:
            base_incar = builder_inputs['parameters']['incar'].copy()
        else:
            base_incar = {}

        # Load code
        code = orm.load_code(code_label)

        # Create stage task
        stage_task = wg.add_task(
            aimd_single_stage_scatter,
            slabs=current_structures,
            temperature=temperature,
            steps=steps,
            code=code,
            base_aimd_parameters=base_incar,
            structure_aimd_overrides=structure_incar_overrides if structure_incar_overrides else None,
            potential_family=builder_inputs.get('potential_family', 'PBE'),
            potential_mapping=builder_inputs.get('potential_mapping', {}),
            options=builder_inputs.get('options', {}),
            kpoints_spacing=builder_inputs.get('kpoints_spacing', 0.5),
            clean_workdir=builder_inputs.get('clean_workdir', False),
            restart_folders=current_remote_folders,
            max_number_jobs=max_concurrent_jobs,
            name=f'stage_{stage_idx}_aimd',
        )
```

**Step 2: Commit override extraction**

```bash
git add teros/core/aimd/workgraph.py
git commit -m "feat(aimd): implement override extraction and priority system

Priority order: matrix > stage > structure > base

Extracts INCAR parameters from all three override levels and
passes merged dict per structure to scatter function."
```

---

## Task 5: Update Docstring for Override Parameters

**Files:**
- Modify: `teros/core/aimd/workgraph.py:67-105`

**Step 1: Update parameter descriptions**

Find the docstring in `build_aimd_workgraph` (starts around line 67). Update the override parameter descriptions:

```python
    """
    Build AIMD workgraph with sequential stages.

    Args:
        structures: {name: StructureData or PK} - input structures
        aimd_stages: [{'temperature': K, 'steps': N}, ...] - sequential stages
        code_label: VASP code label (e.g., 'VASP6.5.0@cluster02')
        builder_inputs: Default builder config for all (structure, stage) combinations
        supercell_specs: {structure_name: [nx, ny, nz]} - optional supercell per structure
        structure_overrides: Per-structure builder overrides.
                           Only 'parameters'/'incar' keys are applied.
                           Other keys (kpoints, options, etc.) are ignored.
                           Format: {structure_name: {'parameters': {'incar': {...}}}}
        stage_overrides: Per-stage builder overrides (0-indexed).
                        Only 'parameters'/'incar' keys are applied.
                        Format: {stage_idx: {'parameters': {'incar': {...}}}}
        matrix_overrides: Per-(structure, stage) builder overrides.
                         Only 'parameters'/'incar' keys are applied.
                         Format: {(structure_name, stage_idx): {'parameters': {'incar': {...}}}}
        max_concurrent_jobs: Limit parallel VASP calculations (None = unlimited)
        name: WorkGraph name

    Returns:
        WorkGraph ready to submit

    Override priority: matrix_overrides > stage_overrides > structure_overrides > builder_inputs

    Note: Only INCAR parameters can be overridden. kpoints_spacing, options,
          potential_mapping, and other builder inputs remain uniform across all structures.

    Example:
        wg = build_aimd_workgraph(
            structures={'slab1': structure1, 'slab2': pk2},
            aimd_stages=[
                {'temperature': 300, 'steps': 100},
                {'temperature': 300, 'steps': 500},
            ],
            code_label='VASP6.5.0@cluster02',
            builder_inputs={
                'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400}},
                'kpoints_spacing': 0.5,
                'potential_family': 'PBE',
                'potential_mapping': {},
                'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
                'clean_workdir': False,
            },
            # Optional overrides
            structure_overrides={
                'slab2': {'parameters': {'incar': {'ENCUT': 500}}}  # slab2 uses ENCUT=500
            },
            stage_overrides={
                1: {'parameters': {'incar': {'PREC': 'Accurate'}}}  # stage 1 uses Accurate
            },
            matrix_overrides={
                ('slab1', 1): {'parameters': {'incar': {'ALGO': 'Fast'}}}  # slab1+stage1 specific
            },
            max_concurrent_jobs=4,
        )
    """
```

**Step 2: Remove old limitation comment**

Find and remove the old "CURRENT LIMITATIONS" section from the docstring (if it exists around lines 88-97).

**Step 3: Commit docstring update**

```bash
git add teros/core/aimd/workgraph.py
git commit -m "docs(aimd): update docstring with functional override system

- Remove 'NOT IMPLEMENTED' warnings
- Document INCAR-only override scope
- Add priority order and example usage
- Clarify what can/cannot be overridden"
```

---

## Task 6: Write Override Tests

**Files:**
- Create: `teros/core/aimd/test_overrides.py`

**Step 1: Write test for structure overrides**

Create new test file:

```python
"""Tests for AIMD override system."""
from aiida import orm, load_profile
from ase.build import bulk
from teros.core.aimd import build_aimd_workgraph

# Load AiiDA profile for tests
load_profile('presto')


def test_structure_overrides():
    """Test per-structure INCAR overrides."""
    # Create test structures
    atoms1 = bulk('Al', 'fcc', a=4.0)
    atoms2 = bulk('Fe', 'bcc', a=2.87)
    struct1 = orm.StructureData(ase=atoms1)
    struct2 = orm.StructureData(ase=atoms2)

    # Base config
    builder_inputs = {
        'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400}},
        'kpoints_spacing': 0.5,
        'potential_family': 'PBE',
        'potential_mapping': {'Al': 'Al', 'Fe': 'Fe'},
        'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        'clean_workdir': False,
    }

    # Override struct1 to use ENCUT=500
    structure_overrides = {
        'struct1': {'parameters': {'incar': {'ENCUT': 500}}}
    }

    # Build workgraph
    wg = build_aimd_workgraph(
        structures={'struct1': struct1, 'struct2': struct2},
        aimd_stages=[{'temperature': 300, 'steps': 10}],
        code_label='VASP6.5.0@cluster02',
        builder_inputs=builder_inputs,
        structure_overrides=structure_overrides,
        name='test_structure_overrides',
    )

    # Verify workgraph was created
    assert wg is not None
    assert wg.name == 'test_structure_overrides'

    # Verify task exists
    assert 'stage_0_aimd' in [task.name for task in wg.tasks]


def test_stage_overrides():
    """Test per-stage INCAR overrides."""
    atoms = bulk('Al', 'fcc', a=4.0)
    struct = orm.StructureData(ase=atoms)

    builder_inputs = {
        'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400}},
        'kpoints_spacing': 0.5,
        'potential_family': 'PBE',
        'potential_mapping': {'Al': 'Al'},
        'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        'clean_workdir': False,
    }

    # Override stage 1 to use PREC=Accurate
    stage_overrides = {
        1: {'parameters': {'incar': {'PREC': 'Accurate'}}}
    }

    wg = build_aimd_workgraph(
        structures={'struct': struct},
        aimd_stages=[
            {'temperature': 300, 'steps': 10},
            {'temperature': 300, 'steps': 20},
        ],
        code_label='VASP6.5.0@cluster02',
        builder_inputs=builder_inputs,
        stage_overrides=stage_overrides,
        name='test_stage_overrides',
    )

    assert wg is not None
    assert 'stage_0_aimd' in [task.name for task in wg.tasks]
    assert 'stage_1_aimd' in [task.name for task in wg.tasks]


def test_matrix_overrides():
    """Test (structure, stage) specific overrides."""
    atoms1 = bulk('Al', 'fcc', a=4.0)
    atoms2 = bulk('Fe', 'bcc', a=2.87)
    struct1 = orm.StructureData(ase=atoms1)
    struct2 = orm.StructureData(ase=atoms2)

    builder_inputs = {
        'parameters': {'incar': {'PREC': 'Normal', 'ENCUT': 400, 'ALGO': 'Normal'}},
        'kpoints_spacing': 0.5,
        'potential_family': 'PBE',
        'potential_mapping': {'Al': 'Al', 'Fe': 'Fe'},
        'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        'clean_workdir': False,
    }

    # Override (struct1, stage 0) specifically
    matrix_overrides = {
        ('struct1', 0): {'parameters': {'incar': {'ALGO': 'Fast'}}}
    }

    wg = build_aimd_workgraph(
        structures={'struct1': struct1, 'struct2': struct2},
        aimd_stages=[{'temperature': 300, 'steps': 10}],
        code_label='VASP6.5.0@cluster02',
        builder_inputs=builder_inputs,
        matrix_overrides=matrix_overrides,
        name='test_matrix_overrides',
    )

    assert wg is not None


def test_override_priority():
    """Test override priority: matrix > stage > structure > base."""
    atoms = bulk('Al', 'fcc', a=4.0)
    struct = orm.StructureData(ase=atoms)

    builder_inputs = {
        'parameters': {'incar': {'ENCUT': 300}},  # base
        'kpoints_spacing': 0.5,
        'potential_family': 'PBE',
        'potential_mapping': {'Al': 'Al'},
        'options': {'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
        'clean_workdir': False,
    }

    # All levels defined - matrix should win
    structure_overrides = {'struct': {'parameters': {'incar': {'ENCUT': 400}}}}
    stage_overrides = {0: {'parameters': {'incar': {'ENCUT': 500}}}}
    matrix_overrides = {('struct', 0): {'parameters': {'incar': {'ENCUT': 600}}}}

    wg = build_aimd_workgraph(
        structures={'struct': struct},
        aimd_stages=[{'temperature': 300, 'steps': 10}],
        code_label='VASP6.5.0@cluster02',
        builder_inputs=builder_inputs,
        structure_overrides=structure_overrides,
        stage_overrides=stage_overrides,
        matrix_overrides=matrix_overrides,
        name='test_priority',
    )

    assert wg is not None
```

**Step 2: Run tests**

```bash
pytest teros/core/aimd/test_overrides.py -v
```

Expected: All 4 tests PASS

**Step 3: Commit tests**

```bash
git add teros/core/aimd/test_overrides.py
git commit -m "test(aimd): add override system tests

Tests:
- test_structure_overrides: per-structure INCAR
- test_stage_overrides: per-stage INCAR
- test_matrix_overrides: specific (structure, stage)
- test_override_priority: verify matrix > stage > structure > base"
```

---

## Task 7: Verify Existing Tests Still Pass

**Files:**
- Test: `teros/core/aimd/test_*.py`

**Step 1: Run all AIMD module tests**

```bash
pytest teros/core/aimd/test_*.py -v
```

Expected: All tests PASS (15 existing + 4 new = 19 total)

**Step 2: If any tests fail**

Check error messages. Most likely issues:
- Import errors (check module __init__.py)
- Signature mismatches (verify all parameter renames)

Fix any failures before proceeding.

---

## Task 8: Update README Documentation

**Files:**
- Modify: `teros/core/aimd/README.md:176-186`

**Step 1: Update "Current Limitations" section**

Find the "Current Limitations" section (around line 176) and replace with:

```markdown
## Override System

The override parameters (`structure_overrides`, `stage_overrides`, `matrix_overrides`)
support **INCAR parameter customization**.

### Supported Parameters

Any VASP INCAR tag can be overridden:
- `ENCUT`, `PREC`, `EDIFF`, `EDIFFG`
- `ALGO`, `ISIF`, `IBRION`
- `NCORE`, `KPAR`, `LREAL`
- Any other INCAR setting

### Not Supported

These builder inputs remain uniform across all structures:
- `kpoints_spacing` - K-points grid density
- `options` - Scheduler settings (num_cores, walltime, etc.)
- `potential_mapping` - Element to pseudopotential mapping
- `potential_family` - Pseudopotential family name

**Workaround:** Use separate `build_aimd_workgraph()` calls for structures needing different kpoints/options.

### Priority Order

When multiple override levels are specified:
1. **Matrix overrides** (highest): `{(structure_name, stage_idx): {...}}`
2. **Stage overrides**: `{stage_idx: {...}}`
3. **Structure overrides**: `{structure_name: {...}}`
4. **Base parameters** (lowest): `builder_inputs['parameters']['incar']`

### Example

```python
wg = build_aimd_workgraph(
    structures={'slab1': s1, 'slab2': s2},
    aimd_stages=[
        {'temperature': 300, 'steps': 100},
        {'temperature': 300, 'steps': 500},
    ],
    code_label='VASP6.5.0@cluster02',
    builder_inputs={
        'parameters': {'incar': {'ENCUT': 400, 'PREC': 'Normal'}},  # Base
        # ... other parameters
    },

    # slab2 needs higher cutoff everywhere
    structure_overrides={
        'slab2': {'parameters': {'incar': {'ENCUT': 500}}}
    },

    # Stage 1 (production) needs tighter convergence
    stage_overrides={
        1: {'parameters': {'incar': {'EDIFF': 1e-7}}}
    },

    # slab1 + stage 1 needs special algorithm
    matrix_overrides={
        ('slab1', 1): {'parameters': {'incar': {'ALGO': 'All'}}}
    },
)
```

**Result:**
- slab1, stage 0: `ENCUT=400, PREC=Normal` (base)
- slab1, stage 1: `ENCUT=400, PREC=Normal, EDIFF=1e-7, ALGO=All` (stage + matrix)
- slab2, stage 0: `ENCUT=500, PREC=Normal` (structure override)
- slab2, stage 1: `ENCUT=500, PREC=Normal, EDIFF=1e-7` (structure + stage)
```

**Step 2: Commit README update**

```bash
git add teros/core/aimd/README.md
git commit -m "docs(aimd): document functional override system in README

- Replace limitation section with override documentation
- Document supported/unsupported parameters
- Add priority order explanation
- Include comprehensive example with results"
```

---

## Task 9: Create Integration Test Example

**Files:**
- Create: `examples/vasp/step_19_aimd_with_overrides.py`

**Step 1: Create example script**

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 19: AIMD with Per-Structure Overrides

Demonstrates the override system in standalone AIMD module:
- Structure-level overrides
- Stage-level overrides
- Matrix-level overrides
- Priority order verification

Material: Ag (silver bulk)
AIMD: 2-stage sequence
Structures: 2 identical structures with different INCAR settings

Usage:
    source ~/envs/aiida/bin/activate
    python step_19_aimd_with_overrides.py
"""

import sys
import os
from aiida import load_profile, orm
from teros.core.aimd import build_aimd_workgraph
from ase.io import read


def main():
    """Step 19: Test AIMD override system."""

    print("\n" + "="*70)
    print("STEP 19: AIMD WITH PER-STRUCTURE OVERRIDES")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='presto')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')
    ag_cif = os.path.join(structures_dir, 'Ag.cif')

    print(f"\n2. Loading structures:")
    print(f"   Bulk: {ag_cif}")

    # Load structures using ASE
    ag_ase = read(ag_cif)
    ag_structure1 = orm.StructureData(ase=ag_ase)
    ag_structure2 = orm.StructureData(ase=ag_ase)

    print(f"   ✓ Loaded 2 Ag structures")
    print(f"     - Structure 1: {len(ag_structure1.sites)} atoms")
    print(f"     - Structure 2: {len(ag_structure2.sites)} atoms")

    # Code configuration
    code_label = 'VASP6.5.0@cluster02'
    potential_family = 'PBE'

    print(f"\n3. VASP configuration:")
    print(f"   Code: {code_label}")
    print(f"   Potential family: {potential_family}")

    # AIMD stages
    print("\n4. AIMD configuration:")
    aimd_stages = [
        {'temperature': 300, 'steps': 50},   # Equilibration
        {'temperature': 300, 'steps': 100},  # Production
    ]

    for i, stage in enumerate(aimd_stages):
        print(f"   Stage {i}: {stage['temperature']} K, {stage['steps']} steps")

    # Base builder inputs
    builder_inputs = {
        'parameters': {
            'incar': {
                # Basic settings (BASE)
                'PREC': 'Normal',
                'ENCUT': 400,
                'EDIFF': 1e-5,
                'ISMEAR': 0,
                'SIGMA': 0.05,
                'ALGO': 'Normal',
                'LREAL': 'Auto',

                # AIMD settings
                'IBRION': 0,
                'MDALGO': 2,
                'POTIM': 2.0,
                'SMASS': 0.0,

                # Output control
                'LWAVE': False,
                'LCHARG': False,
            }
        },
        'kpoints_spacing': 0.5,
        'potential_family': potential_family,
        'potential_mapping': {'Ag': 'Ag'},
        'options': {
            'resources': {
                'num_machines': 1,
                'num_cores_per_machine': 24,
            },
        },
        'clean_workdir': False,
    }

    print("\n5. Override configuration:")
    print("   Base: ENCUT=400, PREC=Normal, ALGO=Normal")

    # Structure override: ag2 uses higher cutoff
    structure_overrides = {
        'ag2': {'parameters': {'incar': {'ENCUT': 500}}}
    }
    print("   Structure override: ag2 uses ENCUT=500")

    # Stage override: stage 1 (production) uses Accurate
    stage_overrides = {
        1: {'parameters': {'incar': {'PREC': 'Accurate', 'EDIFF': 1e-6}}}
    }
    print("   Stage override: stage 1 uses PREC=Accurate, EDIFF=1e-6")

    # Matrix override: (ag1, stage 1) uses Fast algorithm
    matrix_overrides = {
        ('ag1', 1): {'parameters': {'incar': {'ALGO': 'Fast'}}}
    }
    print("   Matrix override: (ag1, stage 1) uses ALGO=Fast")

    print("\n6. Expected INCAR values:")
    print("   ag1, stage 0: ENCUT=400, PREC=Normal, ALGO=Normal")
    print("   ag1, stage 1: ENCUT=400, PREC=Accurate, ALGO=Fast, EDIFF=1e-6")
    print("   ag2, stage 0: ENCUT=500, PREC=Normal, ALGO=Normal")
    print("   ag2, stage 1: ENCUT=500, PREC=Accurate, ALGO=Normal, EDIFF=1e-6")

    print("\n7. Building workgraph with overrides...")

    # Build AIMD workgraph
    wg = build_aimd_workgraph(
        # Input structures
        structures={
            'ag1': ag_structure1,
            'ag2': ag_structure2,
        },

        # AIMD stages
        aimd_stages=aimd_stages,

        # Code configuration
        code_label=code_label,

        # Base builder inputs
        builder_inputs=builder_inputs,

        # Override system (now functional!)
        structure_overrides=structure_overrides,
        stage_overrides=stage_overrides,
        matrix_overrides=matrix_overrides,

        # Concurrency control
        max_concurrent_jobs=2,

        # Workgraph name
        name='Step19_AIMD_Overrides_Ag',
    )

    print("   ✓ WorkGraph built successfully")

    # Submit
    print("\n8. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 19 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"\nMonitor with:")
    print(f"  verdi process show {wg.pk}")
    print(f"\nVerify INCAR overrides:")
    print(f"  # Find VASP calculation PKs")
    print(f"  verdi process show {wg.pk}  # Look for VaspWorkChain PKs")
    print(f"  # Check INCAR for each calculation")
    print(f"  verdi calcjob inputcat <VASP_PK> INCAR | grep -E 'ENCUT|PREC|ALGO|EDIFF'")
    print(f"\nExpected to see different INCAR values per structure/stage!")
    print(f"{'='*70}\n")

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

**Step 2: Make executable and test syntax**

```bash
chmod +x examples/vasp/step_19_aimd_with_overrides.py
python -m py_compile examples/vasp/step_19_aimd_with_overrides.py
```

Expected: No syntax errors

**Step 3: Commit example**

```bash
git add examples/vasp/step_19_aimd_with_overrides.py
git commit -m "feat(examples): add AIMD override system demonstration

Example demonstrates:
- Structure-level INCAR overrides
- Stage-level INCAR overrides
- Matrix-level INCAR overrides
- Priority order verification
- How to verify overrides applied correctly"
```

---

## Task 10: Run Integration Test

**Files:**
- Test: `examples/vasp/step_19_aimd_with_overrides.py`

**Step 1: Clear Python cache and restart daemon**

```bash
find teros/core -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null
find teros/core -name "*.pyc" -delete 2>/dev/null
verdi daemon restart
```

Expected: Daemon restarts successfully

**Step 2: Run example script**

```bash
source ~/envs/aiida/bin/activate
/home/thiagotd/envs/aiida/bin/python examples/vasp/step_19_aimd_with_overrides.py
```

Expected: Script completes, prints WorkGraph PK

**Step 3: Wait and check status**

```bash
sleep 30
verdi process show <PK_from_output>
```

Expected: WorkGraph state should be "Waiting" or "Running", not "Failed"

**Step 4: Verify override application (once calculations start)**

```bash
# Get stage_0 WorkGraph PK from main WorkGraph
verdi process show <main_PK>

# Get VASP calculation PKs from stage_0
verdi process show <stage_0_PK>

# Check INCAR for ag1 and ag2 calculations
verdi calcjob inputcat <VASP_ag1_PK> INCAR | grep -E 'ENCUT|PREC|ALGO'
verdi calcjob inputcat <VASP_ag2_PK> INCAR | grep -E 'ENCUT|PREC|ALGO'
```

Expected for ag1 stage 0: `ENCUT = 400`, `PREC = Normal`, `ALGO = Normal`
Expected for ag2 stage 0: `ENCUT = 500`, `PREC = Normal`, `ALGO = Normal`

---

## Task 11: Final Verification and Documentation

**Files:**
- Verify: All changes committed
- Update: `teros/core/aimd/README.md` if needed

**Step 1: Run all tests one final time**

```bash
pytest teros/core/aimd/test_*.py -v
```

Expected: All 19 tests PASS

**Step 2: Verify git status**

```bash
git status
```

Expected: Clean working directory or only untracked files

**Step 3: Review commit history**

```bash
git log --oneline -11
```

Expected: See all 11 commits from this implementation

**Step 4: Update main README if needed**

Check if `teros/core/aimd/README.md` needs any final clarifications based on testing experience.

---

## Success Criteria Checklist

- [ ] All 19 tests pass (15 existing + 4 new)
- [ ] Example script runs without errors
- [ ] WorkGraph starts successfully (not Failed immediately)
- [ ] INCAR files show different values for ag1 vs ag2
- [ ] Documentation updated (README, docstrings)
- [ ] All changes committed (11 commits total)
- [ ] No Python cache files remaining

## Troubleshooting

### Issue: Tests fail with import error

**Solution:** Check `teros/core/aimd/__init__.py` exports the function correctly

### Issue: WorkGraph fails immediately

**Solution:** Check `verdi process report <PK>` for error message. Most likely:
- Signature mismatch in scatter function call
- Missing parameter in workgraph.py

### Issue: INCAR values not overridden

**Solution:** Check that override extraction logic correctly navigates nested dict structure. Print `structure_incar_overrides` to debug.

### Issue: Example script can't find code

**Solution:** Verify code label exists: `verdi code list | grep VASP6.5.0@cluster02`
