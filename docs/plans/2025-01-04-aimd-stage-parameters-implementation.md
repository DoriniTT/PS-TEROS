# AIMD Stage-Specific Parameters Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace temperature/steps with VASP-native AIMD parameters in stage definitions

**Architecture:** Modify prepare_aimd_parameters() to accept stage_config dict instead of temperature/steps. Update aimd_single_stage_scatter() signature accordingly. Maintain priority: override system > stage params > base INCAR.

**Tech Stack:** Python 3.13, AiiDA, AiiDA-WorkGraph, pytest

**Design:** `docs/plans/2025-01-04-aimd-stage-parameters-design.md`

---

## Task 1: Write Failing Test for New prepare_aimd_parameters() Signature

**Files:**
- Create: `teros/core/test_aimd_functions.py`

**Step 1: Write failing test for required parameters**

Create new test file:

```python
"""Tests for AIMD functions."""
import pytest
from teros.core.aimd_functions import prepare_aimd_parameters


def test_prepare_aimd_parameters_required_only():
    """Test prepare_aimd_parameters with only required TEBEG and NSW."""
    base_params = {
        'PREC': 'Normal',
        'ENCUT': 400,
        'EDIFF': 1e-5,
    }

    stage_config = {
        'TEBEG': 300,
        'NSW': 100,
    }

    result = prepare_aimd_parameters(base_params, stage_config)

    # Should have base parameters
    assert result['PREC'] == 'Normal'
    assert result['ENCUT'] == 400
    assert result['EDIFF'] == 1e-5

    # Should have AIMD parameters from stage
    assert result['TEBEG'] == 300
    assert result['TEEND'] == 300  # Defaults to TEBEG
    assert result['NSW'] == 100
    assert result['IBRION'] == 0  # Forced to 0 for MD


def test_prepare_aimd_parameters_with_optional():
    """Test prepare_aimd_parameters with optional AIMD parameters."""
    base_params = {
        'PREC': 'Normal',
        'POTIM': 1.0,  # Base timestep
    }

    stage_config = {
        'TEBEG': 300,
        'NSW': 100,
        'TEEND': 400,  # Different end temp
        'POTIM': 2.0,  # Override base
        'MDALGO': 2,
        'SMASS': 0.0,
    }

    result = prepare_aimd_parameters(base_params, stage_config)

    # Stage AIMD params should override base
    assert result['TEBEG'] == 300
    assert result['TEEND'] == 400
    assert result['NSW'] == 100
    assert result['POTIM'] == 2.0  # Overridden
    assert result['MDALGO'] == 2
    assert result['SMASS'] == 0.0


def test_prepare_aimd_parameters_missing_tebeg():
    """Test that ValueError is raised if TEBEG missing."""
    base_params = {'PREC': 'Normal'}
    stage_config = {'NSW': 100}

    with pytest.raises(ValueError, match="must contain 'TEBEG' and 'NSW'"):
        prepare_aimd_parameters(base_params, stage_config)


def test_prepare_aimd_parameters_missing_nsw():
    """Test that ValueError is raised if NSW missing."""
    base_params = {'PREC': 'Normal'}
    stage_config = {'TEBEG': 300}

    with pytest.raises(ValueError, match="must contain 'TEBEG' and 'NSW'"):
        prepare_aimd_parameters(base_params, stage_config)
```

**Step 2: Run test to verify it fails**

```bash
pytest teros/core/test_aimd_functions.py -v
```

Expected output:
```
FAILED - TypeError: prepare_aimd_parameters() takes 3 positional arguments but 2 were given
```

**Step 3: Commit failing test**

```bash
git add teros/core/test_aimd_functions.py
git commit -m "test(aimd): add tests for new prepare_aimd_parameters signature

- Test required parameters (TEBEG, NSW)
- Test optional parameters (TEEND, POTIM, MDALGO, SMASS)
- Test validation (missing TEBEG or NSW raises ValueError)
- Tests currently fail (signature not yet updated)"
```

---

## Task 2: Update prepare_aimd_parameters() Implementation

**Files:**
- Modify: `teros/core/aimd_functions.py:21-55`

**Step 1: Update function signature and implementation**

Find `prepare_aimd_parameters()` (around line 21) and replace:

```python
def prepare_aimd_parameters(base_aimd_parameters: dict, stage_config: dict) -> dict:
    """
    Build AIMD INCAR from base parameters and stage-specific AIMD config.

    Args:
        base_aimd_parameters: Base INCAR parameters
        stage_config: Stage configuration dict with required TEBEG, NSW
                     and optional TEEND, POTIM, MDALGO, SMASS parameters

    Returns:
        Complete INCAR dict for this AIMD stage

    Raises:
        ValueError: If TEBEG or NSW not in stage_config
    """
    # Start with base parameters
    aimd_incar = base_aimd_parameters.copy()

    # Validate required parameters
    if 'TEBEG' not in stage_config or 'NSW' not in stage_config:
        raise ValueError("aimd_stages dict must contain 'TEBEG' and 'NSW' parameters")

    # Extract required parameters
    tebeg = stage_config['TEBEG']
    nsw = stage_config['NSW']

    # TEEND defaults to TEBEG (constant temperature MD)
    teend = stage_config.get('TEEND', tebeg)

    # Set temperature and steps
    aimd_incar['TEBEG'] = tebeg
    aimd_incar['TEEND'] = teend
    aimd_incar['NSW'] = nsw

    # Set optional AIMD parameters if provided (override base)
    for param in ['POTIM', 'MDALGO', 'SMASS']:
        if param in stage_config:
            aimd_incar[param] = stage_config[param]

    # Ensure IBRION=0 for MD mode
    aimd_incar['IBRION'] = 0

    return aimd_incar
```

**Step 2: Run tests to verify they pass**

```bash
pytest teros/core/test_aimd_functions.py -v
```

Expected output:
```
test_prepare_aimd_parameters_required_only PASSED
test_prepare_aimd_parameters_with_optional PASSED
test_prepare_aimd_parameters_missing_tebeg PASSED
test_prepare_aimd_parameters_missing_nsw PASSED
======================== 4 passed in 0.5s ========================
```

**Step 3: Commit implementation**

```bash
git add teros/core/aimd_functions.py
git commit -m "feat(aimd): update prepare_aimd_parameters to use stage_config

- Replace temperature/steps parameters with stage_config dict
- Support TEBEG, NSW (required) and TEEND, POTIM, MDALGO, SMASS (optional)
- TEEND defaults to TEBEG for constant-temperature MD
- Validate required parameters, raise ValueError if missing
- All 4 new tests pass"
```

---

## Task 3: Update aimd_single_stage_scatter() Signature

**Files:**
- Modify: `teros/core/aimd_functions.py:122-195`

**Step 1: Update function signature**

Find `aimd_single_stage_scatter()` (line 122) and update signature:

```python
@task.graph
def aimd_single_stage_scatter(
    slabs: t.Annotated[dict[str, orm.StructureData], dynamic(orm.StructureData)],
    stage_config: dict,  # CHANGED: replaces temperature and steps
    code: orm.Code,
    base_aimd_parameters: dict,
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
    restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = {},
    structure_aimd_overrides: dict[str, dict] = None,
    max_number_jobs: int = None,
) -> t.Annotated[dict, namespace(structures=dynamic(orm.StructureData), remote_folders=dynamic(orm.RemoteData), energies=dynamic(orm.Float))]:
```

**Step 2: Update docstring**

Update docstring (lines ~137-166):

```python
    """
    Run single AIMD stage on all slabs in parallel using scatter-gather pattern.

    This function handles ONE AIMD stage for all slabs.
    Call it multiple times sequentially to build multi-stage AIMD workflows.

    Args:
        slabs: Dictionary of slab structures to run AIMD on
        stage_config: Stage AIMD configuration dict with TEBEG, NSW (required)
                     and TEEND, POTIM, MDALGO, SMASS (optional)
        code: VASP code
        base_aimd_parameters: Base AIMD INCAR parameters
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

**Step 3: Update scatter loop to use new signature**

Find the scatter loop (line ~185-194) and update:

```python
    # Scatter: create AIMD task for each slab (runs in parallel)
    for slab_label, slab_structure in slabs.items():
        # Step 1: Apply stage AIMD parameters to base
        stage_params = prepare_aimd_parameters(base_aimd_parameters, stage_config)

        # Step 2: Merge with structure-specific overrides
        if structure_aimd_overrides and slab_label in structure_aimd_overrides:
            # Shallow merge: structure overrides take precedence
            merged_params = {**stage_params, **structure_aimd_overrides[slab_label]}
        else:
            # No override for this structure
            merged_params = stage_params
```

**Step 4: Commit signature update**

```bash
git add teros/core/aimd_functions.py
git commit -m "refactor(aimd): update aimd_single_stage_scatter signature

- Replace temperature and steps with stage_config dict
- Update docstring to document new parameter
- Update scatter loop: apply stage params then override system
- Maintains priority: override system > stage params > base"
```

---

## Task 4: Update build_aimd_workgraph() Stage Loop

**Files:**
- Modify: `teros/core/aimd/workgraph.py:175-235`

**Step 1: Update stage loop**

Find the stage loop in `build_aimd_workgraph()` (line ~175) and replace:

```python
    for stage_idx, stage_config in enumerate(aimd_stages):
        # Validate required AIMD parameters
        if 'TEBEG' not in stage_config or 'NSW' not in stage_config:
            raise ValueError(
                f"Stage {stage_idx}: aimd_stages must contain 'TEBEG' and 'NSW'. "
                f"Got: {list(stage_config.keys())}"
            )

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
            stage_config=stage_config,  # CHANGED: pass entire config
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

**Step 2: Commit standalone module update**

```bash
git add teros/core/aimd/workgraph.py
git commit -m "feat(aimd): update build_aimd_workgraph to use new stage_config format

- Add validation for required TEBEG and NSW parameters
- Pass entire stage_config dict to scatter function
- Clear error message shows available keys if validation fails"
```

---

## Task 5: Update Main Workflow Caller (if applicable)

**Files:**
- Check: `teros/core/workgraph.py:1742-1759`

**Step 1: Search for aimd_single_stage_scatter calls**

```bash
grep -n "aimd_single_stage_scatter" teros/core/workgraph.py
```

**Step 2: If found, update the call**

If the main workflow uses aimd_single_stage_scatter directly (unlikely based on current architecture), update it to use stage_config instead of temperature/steps.

Based on current code review, the main workflow uses the AIMD module indirectly, so no changes needed.

**Step 3: Verify no direct calls exist**

```bash
grep -r "temperature.*steps" teros/core/workgraph.py | grep -i aimd
```

Expected: No matches (main workflow doesn't call scatter directly)

---

## Task 6: Update Example Script step_18

**Files:**
- Modify: `examples/vasp/step_18_aimd_standalone.py:66-69`

**Step 1: Update aimd_stages definition**

Find the aimd_stages definition (line ~66) and replace:

```python
    # AIMD stages
    print("\n4. AIMD configuration:")
    aimd_stages = [
        {'TEBEG': 300, 'NSW': 50},    # Equilibration
        {'TEBEG': 300, 'NSW': 100},   # Production
    ]

    for i, stage in enumerate(aimd_stages):
        print(f"   Stage {i}: TEBEG={stage['TEBEG']} K, NSW={stage['NSW']} steps")
```

**Step 2: Add supercell comment**

Add comment near supercell section (around line 138):

```python
    print("\n7. Building workgraph...")

    # Optional: Create supercells before AIMD
    # Uncomment to create 2x2x1 supercell of the slab:
    # supercell_specs = {'ag_100_slab': [2, 2, 1]}

    # Build AIMD workgraph
    wg = build_aimd_workgraph(
```

**Step 3: Run syntax check**

```bash
python -m py_compile examples/vasp/step_18_aimd_standalone.py
```

Expected: No output (syntax OK)

**Step 4: Commit step_18 update**

```bash
git add examples/vasp/step_18_aimd_standalone.py
git commit -m "refactor(examples): update step_18 to use new AIMD stage format

- Replace temperature/steps with TEBEG/NSW
- Update print statements to show VASP parameter names
- Add supercell_specs comment explaining optional feature"
```

---

## Task 7: Update Example Script step_19

**Files:**
- Modify: `examples/vasp/step_19_aimd_with_overrides.py:66-135`

**Step 1: Update aimd_stages definition**

Find aimd_stages (line ~66) and replace:

```python
    # AIMD stages
    print("\n4. AIMD configuration:")
    aimd_stages = [
        {'TEBEG': 300, 'NSW': 50},    # Equilibration
        {'TEBEG': 300, 'NSW': 100},   # Production
    ]

    for i, stage in enumerate(aimd_stages):
        print(f"   Stage {i}: TEBEG={stage['TEBEG']} K, NSW={stage['NSW']} steps")
```

**Step 2: Update expected INCAR values section**

Find the expected values print section (line ~131) and update:

```python
    print("\n6. Expected INCAR values:")
    print("   ag1, stage 0: ENCUT=400, PREC=Normal, ALGO=Normal, TEBEG=300, NSW=50")
    print("   ag1, stage 1: ENCUT=400, PREC=Accurate, ALGO=Fast, EDIFF=1e-6, TEBEG=300, NSW=100")
    print("   ag2, stage 0: ENCUT=500, PREC=Normal, ALGO=Normal, TEBEG=300, NSW=50")
    print("   ag2, stage 1: ENCUT=500, PREC=Accurate, ALGO=Normal, EDIFF=1e-6, TEBEG=300, NSW=100")
```

**Step 3: Add supercell comment**

Add comment in the workgraph building section:

```python
    print("\n7. Building workgraph with overrides...")

    # Optional: Create supercells before AIMD
    # supercell_specs={'ag1': [2, 2, 1], 'ag2': [2, 2, 1]}

    # Build AIMD workgraph
    wg = build_aimd_workgraph(
```

**Step 4: Run syntax check**

```bash
python -m py_compile examples/vasp/step_19_aimd_with_overrides.py
```

Expected: No output (syntax OK)

**Step 5: Commit step_19 update**

```bash
git add examples/vasp/step_19_aimd_with_overrides.py
git commit -m "refactor(examples): update step_19 to use new AIMD stage format

- Replace temperature/steps with TEBEG/NSW
- Update expected INCAR values to include TEBEG/NSW
- Add supercell_specs comment"
```

---

## Task 8: Update Unit Tests

**Files:**
- Modify: `teros/core/aimd/test_overrides.py:66-142`

**Step 1: Update all aimd_stages definitions in test file**

Find all `aimd_stages` definitions and replace:

```python
# Line ~68 (test_structure_overrides)
aimd_stages=[{'TEBEG': 300, 'NSW': 10}],

# Line ~105 (test_stage_overrides)
aimd_stages=[
    {'TEBEG': 300, 'NSW': 10},
    {'TEBEG': 300, 'NSW': 20},
],

# Line ~138 (test_matrix_overrides)
aimd_stages=[{'TEBEG': 300, 'NSW': 10}],

# Line ~175 (test_override_priority)
aimd_stages=[{'TEBEG': 300, 'NSW': 10}],
```

**Step 2: Run override tests**

```bash
pytest teros/core/aimd/test_overrides.py -v
```

Expected output:
```
test_structure_overrides PASSED
test_stage_overrides PASSED
test_matrix_overrides PASSED
test_override_priority PASSED
======================== 4 passed in 2.5s ========================
```

**Step 3: Commit test updates**

```bash
git add teros/core/aimd/test_overrides.py
git commit -m "test(aimd): update override tests to use new stage format

- Replace temperature/steps with TEBEG/NSW in all tests
- All 4 override tests pass"
```

---

## Task 9: Add Validation Tests

**Files:**
- Modify: `teros/core/aimd/test_overrides.py` (add new tests)

**Step 1: Add test for missing TEBEG**

Append to test file:

```python
def test_build_aimd_workgraph_missing_tebeg():
    """Test that ValueError is raised if TEBEG missing from stage."""
    from aiida import load_profile
    load_profile('presto')

    from ase.build import bulk
    from teros.core.aimd import build_aimd_workgraph

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

    # Missing TEBEG
    aimd_stages = [{'NSW': 100}]

    with pytest.raises(ValueError, match="must contain 'TEBEG' and 'NSW'"):
        build_aimd_workgraph(
            structures={'struct': struct},
            aimd_stages=aimd_stages,
            code_label='VASP-6.5.1@cluster02',
            builder_inputs=builder_inputs,
        )


def test_build_aimd_workgraph_missing_nsw():
    """Test that ValueError is raised if NSW missing from stage."""
    from aiida import load_profile
    load_profile('presto')

    from ase.build import bulk
    from teros.core.aimd import build_aimd_workgraph

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

    # Missing NSW
    aimd_stages = [{'TEBEG': 300}]

    with pytest.raises(ValueError, match="must contain 'TEBEG' and 'NSW'"):
        build_aimd_workgraph(
            structures={'struct': struct},
            aimd_stages=aimd_stages,
            code_label='VASP-6.5.1@cluster02',
            builder_inputs=builder_inputs,
        )
```

**Step 2: Run validation tests**

```bash
pytest teros/core/aimd/test_overrides.py::test_build_aimd_workgraph_missing_tebeg -v
pytest teros/core/aimd/test_overrides.py::test_build_aimd_workgraph_missing_nsw -v
```

Expected: Both tests PASS

**Step 3: Commit validation tests**

```bash
git add teros/core/aimd/test_overrides.py
git commit -m "test(aimd): add validation tests for stage config

- Test missing TEBEG raises ValueError
- Test missing NSW raises ValueError
- Validates error messages are clear"
```

---

## Task 10: Run All AIMD Tests

**Files:**
- Test: `teros/core/aimd/test_*.py`
- Test: `teros/core/test_aimd_functions.py`

**Step 1: Run all AIMD module tests**

```bash
pytest teros/core/aimd/test_*.py teros/core/test_aimd_functions.py -v
```

Expected output:
```
teros/core/test_aimd_functions.py::test_prepare_aimd_parameters_required_only PASSED
teros/core/test_aimd_functions.py::test_prepare_aimd_parameters_with_optional PASSED
teros/core/test_aimd_functions.py::test_prepare_aimd_parameters_missing_tebeg PASSED
teros/core/test_aimd_functions.py::test_prepare_aimd_parameters_missing_nsw PASSED
teros/core/aimd/test_overrides.py::test_structure_overrides PASSED
teros/core/aimd/test_overrides.py::test_stage_overrides PASSED
teros/core/aimd/test_overrides.py::test_matrix_overrides PASSED
teros/core/aimd/test_override_priority PASSED
teros/core/aimd/test_overrides.py::test_build_aimd_workgraph_missing_tebeg PASSED
teros/core/aimd/test_overrides.py::test_build_aimd_workgraph_missing_nsw PASSED
teros/core/aimd/test_tasks.py::test_create_supercell_basic PASSED
teros/core/aimd/test_tasks.py::test_create_supercell_3x3x2 PASSED
teros/core/aimd/test_utils.py::test_validate_stage_sequence_valid PASSED
...
======================== 25 passed in 3.5s ========================
```

**Step 2: If any tests fail**

Check error messages. Most likely issues:
- Import errors (check module paths)
- Signature mismatches (verify all updates complete)

---

## Task 11: Update Module README

**Files:**
- Modify: `teros/core/aimd/README.md:55-85`

**Step 1: Update Quick Start example**

Find the Quick Start section and update aimd_stages:

```markdown
## Quick Start

```python
from teros.core.aimd import build_aimd_workgraph

# Define AIMD stages with VASP-native parameters
aimd_stages = [
    {
        'TEBEG': 300,      # Initial temperature (K) - required
        'NSW': 100,        # Number of MD steps - required
        'POTIM': 2.0,      # Timestep (fs) - optional
    },
    {
        'TEBEG': 300,
        'NSW': 500,
        'POTIM': 1.5,      # Smaller timestep for production
    },
]
```
```

**Step 2: Add new section for stage parameters**

Add after Quick Start section:

```markdown
## AIMD Stage Parameters

Each stage in `aimd_stages` must contain:

**Required**:
- `TEBEG`: Initial temperature (K)
- `NSW`: Number of MD steps

**Optional**:
- `TEEND`: Final temperature (defaults to TEBEG for constant-temp MD)
- `POTIM`: MD timestep in femtoseconds
- `MDALGO`: Thermostat algorithm (1=Andersen, 2=Nosé-Hoover, 3=Langevin)
- `SMASS`: Nosé mass parameter

Missing optional parameters are taken from `builder_inputs['parameters']['incar']` or VASP defaults.

### Examples

**Constant temperature (equilibration + production)**:
```python
aimd_stages = [
    {'TEBEG': 300, 'NSW': 100, 'POTIM': 2.0},   # Equilibration
    {'TEBEG': 300, 'NSW': 500, 'POTIM': 1.5},   # Production
]
```

**Temperature annealing**:
```python
aimd_stages = [
    {'TEBEG': 300, 'TEEND': 500, 'NSW': 200},  # Heat 300K→500K
    {'TEBEG': 500, 'TEEND': 500, 'NSW': 300},  # Hold at 500K
    {'TEBEG': 500, 'TEEND': 300, 'NSW': 200},  # Cool 500K→300K
]
```
```

**Step 3: Update all other examples in README**

Search for `temperature` and `steps` in README and replace with `TEBEG` and `NSW`.

**Step 4: Commit README update**

```bash
git add teros/core/aimd/README.md
git commit -m "docs(aimd): update README for new stage parameter format

- Replace temperature/steps with TEBEG/NSW in examples
- Add new section documenting stage parameters
- Add examples for constant-temp and annealing workflows"
```

---

## Task 12: Update Standalone Module Documentation

**Files:**
- Modify: `docs/aimd_standalone_module.md:66-145`

**Step 1: Update Quick Start example**

Find Quick Start section (line ~66) and update:

```markdown
## Quick Start

### Basic Example: Single Structure, Two Stages

```python
# Define AIMD stages
aimd_stages = [
    {'TEBEG': 300, 'NSW': 100},   # Equilibration
    {'TEBEG': 300, 'NSW': 500},   # Production
]
```
```

**Step 2: Update all examples throughout documentation**

Search and replace in entire file:
- `'temperature': 300` → `'TEBEG': 300`
- `'steps': 100` → `'NSW': 100`

**Step 3: Update Use Case examples**

Update all use case sections (lines ~200-280):

```markdown
### Use Case 1: Temperature Series

```python
for temp in [300, 400, 500, 600]:
    wg = build_aimd_workgraph(
        structures={'slab': structure},
        aimd_stages=[{'TEBEG': temp, 'NSW': 200}],
        # ...
    )
```

### Use Case 3: Equilibration + Production with Different Timesteps

```python
aimd_stages = [
    {'TEBEG': 300, 'NSW': 100, 'POTIM': 2.0},   # Equilibration
    {'TEBEG': 300, 'NSW': 500, 'POTIM': 1.5},   # Production
]
```
```

**Step 4: Add migration note at top**

Add after Overview section:

```markdown
## Breaking Change Notice (2025-01-04)

**IMPORTANT**: The stage definition format changed from `temperature`/`steps` to VASP-native parameter names.

**Before**:
```python
aimd_stages = [{'temperature': 300, 'steps': 100}]
```

**After**:
```python
aimd_stages = [{'TEBEG': 300, 'NSW': 100}]
```

See design document: `docs/plans/2025-01-04-aimd-stage-parameters-design.md`
```

**Step 5: Commit documentation update**

```bash
git add docs/aimd_standalone_module.md
git commit -m "docs(aimd): update standalone module docs for new stage format

- Replace temperature/steps with TEBEG/NSW throughout
- Update all examples and use cases
- Add breaking change migration notice
- Add examples using optional parameters (POTIM, TEEND)"
```

---

## Task 13: Clear Cache and Test Integration

**Files:**
- Test: `examples/vasp/step_19_aimd_with_overrides.py`

**Step 1: Clear Python cache**

```bash
find teros/core -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null
find teros/core -name "*.pyc" -delete 2>/dev/null
```

**Step 2: Restart AiiDA daemon**

```bash
verdi daemon restart
```

Expected: Daemon restarts successfully

**Step 3: Run step_19 integration test**

```bash
source ~/envs/aiida/bin/activate
/home/thiagotd/envs/aiida/bin/python examples/vasp/step_19_aimd_with_overrides.py
```

Expected output:
```
STEP 19: AIMD WITH PER-STRUCTURE OVERRIDES
...
WorkGraph PK: <PK>
```

**Step 4: Wait and verify**

```bash
sleep 30
verdi process show <PK>
```

Expected: WorkGraph state is "Waiting" or "Running", not "Failed"

**Step 5: Verify INCAR contains TEBEG and NSW**

```bash
# Get VASP calculation PK from workgraph
verdi process show <stage_0_PK>

# Check INCAR
verdi calcjob inputcat <VASP_PK> INCAR | grep -E 'TEBEG|NSW|TEEND'
```

Expected output:
```
TEBEG = 300
TEEND = 300
NSW = 50
```

---

## Task 14: Final Verification

**Files:**
- Verify: All changes committed
- Test: All tests passing

**Step 1: Run all tests one final time**

```bash
pytest teros/core/aimd/test_*.py teros/core/test_aimd_functions.py -v
```

Expected: All 25+ tests PASS

**Step 2: Verify git status**

```bash
git status
```

Expected: Clean working directory or only untracked files

**Step 3: Review commit history**

```bash
git log --oneline -14
```

Expected: See all 14 commits from this implementation

---

## Success Criteria Checklist

- [ ] All 25+ tests pass
- [ ] Step 19 integration test runs without errors
- [ ] WorkGraph starts successfully (not Failed immediately)
- [ ] INCAR files contain TEBEG, NSW, TEEND instead of old format
- [ ] Documentation updated (README, standalone module docs)
- [ ] All commits clean and descriptive
- [ ] No Python cache files remaining

## Troubleshooting

### Issue: Tests fail with "aimd_stages must contain 'TEBEG'"

**Solution**: Update all aimd_stages definitions in tests to use new format

### Issue: Step 19 fails immediately

**Solution**: Check `verdi process report <PK>` for error. Most likely:
- Stage validation error (missing TEBEG or NSW)
- Check aimd_stages definition in script

### Issue: INCAR still shows old temperature/steps behavior

**Solution**:
- Verify prepare_aimd_parameters() signature updated
- Clear Python cache and restart daemon
- Check that stage_config is being passed correctly

---

## Execution Notes

This implementation:
- **Breaking change**: All existing code must update aimd_stages format
- **Clear errors**: Validation provides helpful error messages
- **TDD**: All tests written before or alongside implementation
- **Documentation**: Comprehensive updates to all docs
- **Integration tested**: Step 19 verifies end-to-end functionality
