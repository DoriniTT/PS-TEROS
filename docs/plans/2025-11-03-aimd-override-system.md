# AIMD Override System Implementation

**Date:** 2025-11-03
**Status:** Design Complete
**Purpose:** Enable per-structure INCAR parameter overrides in standalone AIMD module

## Problem Statement

The standalone AIMD module (`teros/core/aimd/`) accepts override parameters but cannot apply them. All structures in all stages use identical INCAR parameters from `builder_inputs`.

**Root cause:** `aimd_single_stage_scatter()` accepts a single `aimd_parameters` dict and applies it uniformly to all structures in the scatter loop.

**User impact:** Cannot customize ENCUT, PREC, ALGO, or other INCAR settings per structure or stage.

## Solution: Two-Parameter Approach

Modify `aimd_single_stage_scatter()` to accept:
1. `base_aimd_parameters` (dict) - Applied to all structures by default
2. `structure_aimd_overrides` (dict[str, dict]) - Per-structure INCAR replacements

Structures missing from override dict fall back to base parameters.

**Scope:** INCAR parameters only. Other builder inputs (kpoints, options, potential_mapping) remain uniform across structures.

**Compatibility:** Breaking change. Requires updating existing callers.

## Design

### 1. Function Signature Change

**File:** `teros/core/aimd_functions.py`

```python
@task.graph
def aimd_single_stage_scatter(
    slabs: dict[str, orm.StructureData],
    temperature: float,
    steps: int,
    code: orm.Code,
    base_aimd_parameters: dict,  # RENAMED from aimd_parameters
    potential_family: str,
    potential_mapping: dict,
    options: dict,
    kpoints_spacing: float,
    clean_workdir: bool,
    restart_folders: dict[str, orm.RemoteData] = {},
    structure_aimd_overrides: dict[str, dict] = None,  # NEW
    max_number_jobs: int = None,
) -> dict:
    """
    Run single AIMD stage on all slabs in parallel.

    Args:
        base_aimd_parameters: Base AIMD INCAR parameters applied to all structures
        structure_aimd_overrides: Optional per-structure INCAR overrides.
                                 Format: {structure_name: {INCAR_key: value}}
                                 Missing structures use base_aimd_parameters.
    ...
    """
```

### 2. Implementation Logic

**Location:** Inside `aimd_single_stage_scatter()`, scatter loop

```python
for slab_label, slab_structure in slabs.items():
    # Merge base parameters with structure-specific overrides
    if structure_aimd_overrides and slab_label in structure_aimd_overrides:
        merged_params = {**base_aimd_parameters, **structure_aimd_overrides[slab_label]}
    else:
        merged_params = base_aimd_parameters

    # Prepare parameters for this stage (adds TEBEG, TEEND, NSW)
    stage_params = prepare_aimd_parameters(merged_params, temperature, steps)

    # Build VASP inputs (unchanged)
    vasp_inputs = {
        'structure': slab_structure,
        'code': code,
        'parameters': {'incar': stage_params},
        # ... rest unchanged
    }

    aimd_task = VaspTask(**vasp_inputs)
    # ...
```

**Merge behavior:**
- Shallow merge sufficient (INCAR values are primitives)
- Override values completely replace base values for same key
- No deep merging needed

### 3. Standalone Module Integration

**File:** `teros/core/aimd/workgraph.py`

Extract INCAR parameters from the three override levels:

```python
for stage_idx, stage_config in enumerate(aimd_stages):
    temperature = stage_config['temperature']
    steps = stage_config['steps']

    # Build per-structure INCAR overrides for this stage
    structure_incar_overrides = {}

    for struct_name in prepared_structures:
        override = {}

        # Apply structure-level override
        if structure_overrides and struct_name in structure_overrides:
            struct_override = structure_overrides[struct_name]
            if 'parameters' in struct_override and 'incar' in struct_override['parameters']:
                override.update(struct_override['parameters']['incar'])

        # Apply stage-level override
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

        if override:
            structure_incar_overrides[struct_name] = override

    # Extract base INCAR from builder_inputs
    base_incar = builder_inputs.get('parameters', {}).get('incar', {})

    # Create AIMD task
    stage_task = wg.add_task(
        aimd_single_stage_scatter,
        slabs=current_structures,
        temperature=temperature,
        steps=steps,
        code=code,
        base_aimd_parameters=base_incar,
        structure_aimd_overrides=structure_incar_overrides or None,
        # ... other parameters
    )
```

**Priority order (highest to lowest):**
1. Matrix overrides: `(structure_name, stage_idx)`
2. Stage overrides: `stage_idx`
3. Structure overrides: `structure_name`
4. Base: `builder_inputs['parameters']['incar']`

### 4. Documentation Updates

**File:** `teros/core/aimd/workgraph.py` docstring

Remove "NOT IMPLEMENTED" warnings. Update parameter descriptions:

```python
"""
Args:
    structure_overrides: Override builder_inputs per structure.
                        Only 'parameters'/'incar' keys are applied.
                        Other keys ignored (kpoints, options, etc.).
    stage_overrides: Override builder_inputs per stage (0-indexed).
                     Only 'parameters'/'incar' keys are applied.
    matrix_overrides: Override for specific (structure, stage) combinations.
                      Format: {(structure_name, stage_idx): builder_inputs}
                      Only 'parameters'/'incar' keys are applied.

Override priority: matrix_overrides > stage_overrides > structure_overrides > builder_inputs
"""
```

**File:** `teros/core/aimd/README.md`

Update "Current Limitations" section:

```markdown
## Override System

The override parameters (`structure_overrides`, `stage_overrides`, `matrix_overrides`)
support **INCAR parameter customization only**.

**Supported:** ENCUT, PREC, EDIFF, ALGO, ISIF, etc. (any VASP INCAR tag)

**Not supported:** kpoints_spacing, options, potential_mapping, potential_family
(these remain uniform across all structures)

**Priority:** matrix > stage > structure > base
```

## Testing Strategy

### New tests

**File:** `teros/core/aimd/test_workgraph.py` (new)

```python
def test_structure_overrides():
    """Verify per-structure INCAR overrides apply correctly."""
    wg = build_aimd_workgraph(
        structures={'s1': struct1, 's2': struct2},
        aimd_stages=[{'temperature': 300, 'steps': 10}],
        code_label='VASP6.5.0@cluster02',
        builder_inputs={'parameters': {'incar': {'ENCUT': 400}}},
        structure_overrides={'s1': {'parameters': {'incar': {'ENCUT': 500}}}},
    )
    # Verify s1 uses ENCUT=500, s2 uses ENCUT=400

def test_stage_overrides():
    """Verify per-stage INCAR overrides apply correctly."""
    # Stage 0: PREC=Normal, Stage 1: PREC=Accurate

def test_matrix_overrides():
    """Verify (structure, stage) specific overrides apply correctly."""
    # (s1, 0): ALGO=Fast, others: ALGO=Normal

def test_override_priority():
    """Verify priority: matrix > stage > structure > base."""
    # All four levels defined, matrix wins
```

### Existing tests to update

No changes needed. Unit tests for `validate_*` and `merge_*` functions remain valid.

## Migration Impact

### Files requiring updates

1. **`teros/core/workgraph.py`** - Main workflow
   - Line ~1660: Rename `aimd_parameters=` to `base_aimd_parameters=`
   - Add `structure_aimd_overrides=None`

2. **`examples/vasp/step_07_aimd_simulation.py`**
   - Uses main workflow, breaks indirectly
   - Fix by updating workgraph.py

3. **`teros/core/aimd/workgraph.py`**
   - Implement override extraction logic
   - Update parameter passing to scatter function

### Migration guide

```python
# OLD (breaks):
aimd_single_stage_scatter(
    slabs=structures,
    aimd_parameters={'PREC': 'Normal'},
    ...
)

# NEW (required):
aimd_single_stage_scatter(
    slabs=structures,
    base_aimd_parameters={'PREC': 'Normal'},
    structure_aimd_overrides=None,  # or omit
    ...
)
```

## Success Criteria

Implementation succeeds when:

1. All 15 existing unit tests pass
2. New override tests pass (4 new tests)
3. Example workflow runs with overrides applied correctly
4. `verdi calcjob inputcat <PK> INCAR` shows correct ENCUT/PREC per structure
5. Documentation updated (README, docstrings)
6. Step 07 example works after migration

## Implementation Steps

1. Update `aimd_single_stage_scatter()` signature and logic
2. Update `teros/core/workgraph.py` caller
3. Implement override extraction in `teros/core/aimd/workgraph.py`
4. Write new tests
5. Update documentation
6. Verify with running workflow
7. Commit changes

## Non-Goals

- Supporting non-INCAR overrides (kpoints, options, etc.)
- Backward compatibility with old signature
- Automatic parameter validation (user responsible for valid INCAR)
- Per-task resource overrides (use separate workflows instead)
