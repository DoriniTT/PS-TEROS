# AIMD Stage-Specific Parameters Design

**Date**: 2025-01-04
**Status**: Approved
**Goal**: Allow AIMD-specific VASP parameters directly in stage definitions

---

## Problem Statement

The AIMD standalone module currently requires users to specify only `temperature` and `steps` in each stage definition:

```python
aimd_stages = [
    {'temperature': 300, 'steps': 100},
    {'temperature': 400, 'steps': 200},
]
```

This design has limitations:

1. **Temperature range not supported**: Users cannot specify different start/end temperatures (TEBEG ≠ TEEND)
2. **Timestep locked globally**: POTIM must be set in base INCAR, cannot vary between stages
3. **Thermostat settings fixed**: MDALGO and SMASS cannot change between equilibration and production
4. **Not VASP-native**: Uses `temperature` and `steps` instead of actual VASP parameter names

Users must work around these by using the override system, which is verbose for stage-specific AIMD parameters:

```python
stage_overrides = {
    0: {'parameters': {'incar': {'POTIM': 2.0}}},  # Equilibration: larger timestep
    1: {'parameters': {'incar': {'POTIM': 1.5}}},  # Production: smaller timestep
}
```

## Proposed Solution

Replace `temperature` and `steps` with explicit VASP parameter names in stage definitions.

### New Stage Definition Format

```python
aimd_stages = [
    {
        'TEBEG': 300,      # Required: Initial temperature (K)
        'NSW': 100,        # Required: Number of MD steps
        'TEEND': 300,      # Optional: Final temperature (defaults to TEBEG)
        'POTIM': 2.0,      # Optional: Timestep (fs)
        'MDALGO': 2,       # Optional: Thermostat algorithm
        'SMASS': 0.0,      # Optional: Nosé mass
    },
    {
        'TEBEG': 300,
        'NSW': 500,
        'POTIM': 1.5,      # Different timestep for production
    },
]
```

### Required vs Optional Parameters

**Required in every stage**:
- `TEBEG` - Initial temperature (K)
- `NSW` - Number of MD steps

**Optional parameters**:
- `TEEND` - Final temperature (defaults to TEBEG for constant-temperature MD)
- `POTIM` - MD timestep in femtoseconds
- `MDALGO` - Thermostat algorithm (1=Andersen, 2=Nosé-Hoover, 3=Langevin)
- `SMASS` - Nosé mass parameter

Missing optional parameters come from `builder_inputs['parameters']['incar']` or VASP defaults.

### Benefits

1. **VASP-native**: Uses actual VASP parameter names users know from manual
2. **Stage-specific control**: Different timesteps/thermostats per stage without override system
3. **Temperature annealing**: Can specify TEBEG ≠ TEEND for gradual heating/cooling
4. **Clearer code**: `{'TEBEG': 300, 'NSW': 100}` is more explicit than `{'temperature': 300, 'steps': 100}`

---

## Implementation Design

### 1. Function Signature Changes

**File**: `teros/core/aimd_functions.py`

**Current**:
```python
def prepare_aimd_parameters(aimd_parameters: dict, temperature: float, steps: int) -> dict:
```

**New**:
```python
def prepare_aimd_parameters(base_aimd_parameters: dict, stage_config: dict) -> dict:
    """
    Build AIMD INCAR from base parameters and stage-specific AIMD config.

    Args:
        base_aimd_parameters: Base INCAR parameters
        stage_config: Stage configuration dict with TEBEG, NSW (required)
                     and TEEND, POTIM, MDALGO, SMASS (optional)

    Returns:
        Complete INCAR dict for this AIMD stage
    """
```

### 2. Implementation Logic

**Function body**:
```python
def prepare_aimd_parameters(base_aimd_parameters: dict, stage_config: dict) -> dict:
    # Start with base parameters
    aimd_incar = base_aimd_parameters.copy()

    # Validate required parameters
    if 'TEBEG' not in stage_config or 'NSW' not in stage_config:
        raise ValueError("aimd_stages dict must contain 'TEBEG' and 'NSW' parameters")

    # Extract required parameters
    tebeg = stage_config['TEBEG']
    nsw = stage_config['NSW']

    # TEEND defaults to TEBEG (constant temperature)
    teend = stage_config.get('TEEND', tebeg)

    # Set temperature and steps
    aimd_incar['TEBEG'] = tebeg
    aimd_incar['TEEND'] = teend
    aimd_incar['NSW'] = nsw

    # Set optional AIMD parameters if provided
    for param in ['POTIM', 'MDALGO', 'SMASS']:
        if param in stage_config:
            aimd_incar[param] = stage_config[param]

    # Ensure IBRION=0 for MD mode
    aimd_incar['IBRION'] = 0

    return aimd_incar
```

### 3. Integration with Override System

**Priority order** (highest to lowest):
1. Override system (matrix_overrides > stage_overrides > structure_overrides)
2. Stage AIMD parameters (from stage_config)
3. Base INCAR parameters

**In aimd_single_stage_scatter()**:
```python
for slab_label, slab_structure in slabs.items():
    # Step 1: Apply stage AIMD parameters to base
    stage_params = prepare_aimd_parameters(base_aimd_parameters, stage_config)

    # Step 2: Merge with structure-specific overrides
    if structure_aimd_overrides and slab_label in structure_aimd_overrides:
        merged_params = {**stage_params, **structure_aimd_overrides[slab_label]}
    else:
        merged_params = stage_params

    # merged_params now contains final INCAR values
```

**Example priority resolution**:
```python
# Base INCAR: POTIM=1.5
# stage_config: POTIM=2.0
# matrix_override: POTIM=3.0
# Result: POTIM=3.0 (override wins)

# Base INCAR: POTIM=1.5
# stage_config: POTIM=2.0
# No override
# Result: POTIM=2.0 (stage wins)

# Base INCAR: POTIM=1.5
# stage_config: (no POTIM)
# No override
# Result: POTIM=1.5 (base wins)
```

### 4. Scatter Function Signature Change

**File**: `teros/core/aimd_functions.py`

**Current**:
```python
def aimd_single_stage_scatter(
    slabs: dict[str, orm.StructureData],
    temperature: float,
    steps: int,
    code: orm.Code,
    base_aimd_parameters: dict,
    # ... other params
)
```

**New**:
```python
def aimd_single_stage_scatter(
    slabs: dict[str, orm.StructureData],
    stage_config: dict,  # Replaces temperature and steps
    code: orm.Code,
    base_aimd_parameters: dict,
    # ... other params unchanged
)
```

### 5. Standalone Module Changes

**File**: `teros/core/aimd/workgraph.py`

**Current stage loop**:
```python
for stage_idx, stage_config in enumerate(aimd_stages):
    temperature = stage_config['temperature']
    steps = stage_config['steps']

    stage_task = wg.add_task(
        aimd_single_stage_scatter,
        temperature=temperature,
        steps=steps,
        # ...
    )
```

**New stage loop**:
```python
for stage_idx, stage_config in enumerate(aimd_stages):
    # Validate required parameters
    if 'TEBEG' not in stage_config or 'NSW' not in stage_config:
        raise ValueError(
            f"Stage {stage_idx}: aimd_stages must contain 'TEBEG' and 'NSW'. "
            f"Got: {list(stage_config.keys())}"
        )

    # ... extract overrides (unchanged) ...

    stage_task = wg.add_task(
        aimd_single_stage_scatter,
        stage_config=stage_config,  # Pass entire dict
        # ... other params unchanged
    )
```

---

## Breaking Changes and Migration

### This is a Breaking Change

All existing code using `temperature` and `steps` will break with a clear error message.

### Files Requiring Updates

1. **Example scripts**:
   - `examples/vasp/step_18_aimd_standalone.py`
   - `examples/vasp/step_19_aimd_with_overrides.py`
   - Any user scripts using standalone AIMD module

2. **Unit tests**:
   - `teros/core/aimd/test_overrides.py`
   - `teros/core/aimd/test_tasks.py`
   - Add new tests for AIMD parameter validation

3. **Documentation**:
   - `teros/core/aimd/README.md`
   - `docs/aimd_standalone_module.md`
   - Add migration note

### Migration Examples

**Before**:
```python
aimd_stages = [
    {'temperature': 300, 'steps': 100},
    {'temperature': 400, 'steps': 200},
]
```

**After (minimal)**:
```python
aimd_stages = [
    {'TEBEG': 300, 'NSW': 100},
    {'TEBEG': 400, 'NSW': 200},
]
```

**After (with optional parameters)**:
```python
aimd_stages = [
    {
        'TEBEG': 300,
        'NSW': 100,
        'POTIM': 2.0,     # Larger timestep for equilibration
    },
    {
        'TEBEG': 300,
        'NSW': 500,
        'POTIM': 1.5,     # Smaller timestep for production
        'MDALGO': 2,      # Nosé-Hoover thermostat
    },
]
```

**Temperature annealing example**:
```python
aimd_stages = [
    {'TEBEG': 300, 'TEEND': 500, 'NSW': 200},  # Heat from 300K to 500K
    {'TEBEG': 500, 'TEEND': 500, 'NSW': 300},  # Hold at 500K
    {'TEBEG': 500, 'TEEND': 300, 'NSW': 200},  # Cool to 300K
]
```

### Error Messages

**Missing TEBEG**:
```
ValueError: Stage 0: aimd_stages must contain 'TEBEG' and 'NSW'. Got: ['temperature', 'steps']
```

**Missing NSW**:
```
ValueError: Stage 1: aimd_stages must contain 'TEBEG' and 'NSW'. Got: ['TEBEG']
```

---

## Testing Strategy

### Unit Tests

**New tests in `teros/core/aimd/test_aimd_parameters.py`**:

1. `test_prepare_aimd_parameters_required_only()`
   - Verify TEBEG and NSW are correctly set
   - Verify TEEND defaults to TEBEG

2. `test_prepare_aimd_parameters_with_optional()`
   - Verify POTIM, MDALGO, SMASS are set when provided
   - Verify base INCAR values are preserved when not overridden

3. `test_prepare_aimd_parameters_missing_tebeg()`
   - Verify ValueError raised

4. `test_prepare_aimd_parameters_missing_nsw()`
   - Verify ValueError raised

5. `test_stage_config_validation()`
   - Test validation in build_aimd_workgraph()

### Integration Tests

**Update `examples/vasp/step_19_aimd_with_overrides.py`**:
```python
aimd_stages = [
    {
        'TEBEG': 300,
        'NSW': 50,
        'POTIM': 2.0,
    },
    {
        'TEBEG': 300,
        'NSW': 100,
        'POTIM': 1.5,
    },
]
```

Verify INCAR values:
- Stage 0: POTIM=2.0, TEBEG=300, TEEND=300, NSW=50
- Stage 1: POTIM=1.5, TEBEG=300, TEEND=300, NSW=100

---

## Documentation Updates

### 1. Module README

**File**: `teros/core/aimd/README.md`

Update all examples to use new format. Add section:

```markdown
## AIMD Stage Parameters

Each stage in `aimd_stages` must contain:

**Required**:
- `TEBEG`: Initial temperature (K)
- `NSW`: Number of MD steps

**Optional**:
- `TEEND`: Final temperature (defaults to TEBEG)
- `POTIM`: MD timestep (fs)
- `MDALGO`: Thermostat (1=Andersen, 2=Nosé-Hoover, 3=Langevin)
- `SMASS`: Nosé mass parameter

Other INCAR parameters go in `builder_inputs['parameters']['incar']`.
```

### 2. Standalone Module Documentation

**File**: `docs/aimd_standalone_module.md`

Add migration note at top and update all examples throughout.

### 3. Add Supercell Comment

**File**: `examples/vasp/step_19_aimd_with_overrides.py`

Add comment explaining supercell feature:

```python
# Optional: Create supercells before AIMD
# supercell_specs={'ag1': [2, 2, 1], 'ag2': [3, 3, 1]}
```

---

## Implementation Checklist

- [ ] Modify `prepare_aimd_parameters()` signature and implementation
- [ ] Modify `aimd_single_stage_scatter()` signature
- [ ] Update stage loop in `build_aimd_workgraph()`
- [ ] Update all example scripts
- [ ] Update all unit tests
- [ ] Add new validation tests
- [ ] Update module README
- [ ] Update standalone module documentation
- [ ] Add supercell comment to examples
- [ ] Verify all tests pass
- [ ] Run integration test (step_19)

---

## Alternatives Considered

### Alternative 1: Keep temperature/steps, add optional AIMD params

```python
aimd_stages = [
    {'temperature': 300, 'steps': 100, 'potim': 2.0},
]
```

**Rejected because**:
- Mixes abstraction levels (temperature vs TEBEG)
- Creates ambiguity (what if both temperature and TEBEG provided?)
- Not truly VASP-native

### Alternative 2: Override system wins over stage params

Priority: override > stage > base

**Rejected because**:
- Makes stage params too weak
- Defeats purpose of having stage-specific params
- Stage params should be stronger than base, weaker than overrides

---

## Summary

Replace `temperature` and `steps` with VASP-native parameter names (TEBEG, NSW) in stage definitions. Allow optional AIMD parameters (TEEND, POTIM, MDALGO, SMASS) per stage. Maintain priority: override system > stage params > base INCAR.

Breaking change requiring updates to all examples and tests. Clear error messages guide migration.
