# Input Validation Stage Implementation

## Summary

Successfully extracted the input validation logic (lines 1133-1168 of `build_core_workgraph()`) into a dedicated stage module.

## Files Created

### `/home/thiagotd/git/PS-TEROS/teros/core/stages/input_validation.py`

New module implementing Stage 2 of the workflow construction pipeline.

**Functions:**

1. **`validate_workflow_inputs()`** - Main validation function that:
   - Determines if bulk structure is needed based on workflow flags
   - Validates `structures_dir` and `bulk_name` when required
   - Calls `log_workflow_configuration()` to log the setup
   - Returns a boolean indicating whether bulk is needed

2. **`log_workflow_configuration()`** - Logging function that:
   - Outputs a formatted summary of the workflow preset
   - Shows which features are enabled/disabled with visual indicators

## Files Modified

### `/home/thiagotd/git/PS-TEROS/teros/core/stages/__init__.py`

- Replaced placeholder functions with imports from `input_validation.py`
- Updated `__all__` to export `validate_workflow_inputs` (renamed from `validate_required_inputs` for clarity)
- Updated module docstring to reflect the actual function name

## Function Signature

```python
def validate_workflow_inputs(
    resolved_preset_name: str,
    resolved_flags: Dict[str, bool],
    structures_dir: Optional[str],
    bulk_name: Optional[str],
    metal_name: Optional[str],
    oxygen_name: Optional[str],
    miller_indices: Optional[List[tuple]],
    input_slabs: Optional[Dict[str, Any]],
) -> bool:
    """Validate required inputs and determine if bulk structure is needed."""
```

## Logic Preserved

The extracted code preserves the exact original logic:

1. **Bulk needs determination:**
   ```python
   compute_formation_enthalpy = (metal_name is not None and oxygen_name is not None)
   needs_bulk = (
       compute_formation_enthalpy or
       compute_thermodynamics or
       compute_electronic_properties_bulk or
       (miller_indices is not None and input_slabs is None)
   )
   ```

2. **Validation with detailed error messages:**
   - Lists which features require bulk structure
   - Provides actionable guidance (use `input_slabs` with `aimd_only` preset)

3. **Configuration logging:**
   - Uses checkmark/cross Unicode characters for visual clarity
   - Logs preset name and all resolved flags

## Usage Example

```python
from teros.core.stages import validate_workflow_inputs, StageContext

# After Stage 1 (preset resolution)
needs_bulk = validate_workflow_inputs(
    resolved_preset_name='surface_thermodynamics',
    resolved_flags={
        'relax_slabs': True,
        'compute_thermodynamics': True,
        'compute_cleavage': False,
        'compute_electronic_properties_bulk': False,
    },
    structures_dir='/path/to/structures',
    bulk_name='Ag2O_bulk.vasp',
    metal_name='Ag_bulk.vasp',
    oxygen_name='O2.vasp',
    miller_indices=[(1, 1, 0), (1, 0, 0)],
    input_slabs=None,
)
# Returns: True (bulk is needed)
# Logs workflow configuration
```

## Integration Point

In `build_core_workgraph()`, the original code block:

```python
# ========================================================================
# VALIDATE REQUIRED INPUTS
# ========================================================================

# Determine if bulk structure is needed
compute_formation_enthalpy = (metal_name is not None and oxygen_name is not None)
needs_bulk = (
    compute_formation_enthalpy or
    compute_thermodynamics or
    ...
)
# ... validation and logging ...
```

Can now be replaced with:

```python
from teros.core.stages import validate_workflow_inputs

# Stage 2: Validate inputs
needs_bulk = validate_workflow_inputs(
    resolved_preset_name=resolved_preset_name,
    resolved_flags=resolved_flags,
    structures_dir=structures_dir,
    bulk_name=bulk_name,
    metal_name=metal_name,
    oxygen_name=oxygen_name,
    miller_indices=miller_indices,
    input_slabs=input_slabs,
)
```

## Testing

The module can be unit tested in isolation:

```python
import pytest
from teros.core.stages.input_validation import validate_workflow_inputs

def test_validate_bulk_needed_thermodynamics():
    """Bulk is needed when compute_thermodynamics is True."""
    needs_bulk = validate_workflow_inputs(
        resolved_preset_name='surface_thermodynamics',
        resolved_flags={'compute_thermodynamics': True,
                        'compute_electronic_properties_bulk': False},
        structures_dir='/path',
        bulk_name='bulk.vasp',
        metal_name=None,
        oxygen_name=None,
        miller_indices=None,
        input_slabs={'slab_0': 'mock'},
    )
    assert needs_bulk is True

def test_validate_bulk_not_needed_aimd_only():
    """Bulk is not needed for aimd_only with input_slabs."""
    needs_bulk = validate_workflow_inputs(
        resolved_preset_name='aimd_only',
        resolved_flags={'compute_thermodynamics': False,
                        'compute_electronic_properties_bulk': False},
        structures_dir=None,
        bulk_name=None,
        metal_name=None,
        oxygen_name=None,
        miller_indices=None,
        input_slabs={'slab_0': 'mock'},
    )
    assert needs_bulk is False

def test_validate_missing_bulk_raises():
    """Missing bulk inputs raises ValueError with details."""
    with pytest.raises(ValueError) as exc_info:
        validate_workflow_inputs(
            resolved_preset_name='surface_thermodynamics',
            resolved_flags={'compute_thermodynamics': True,
                            'compute_electronic_properties_bulk': False},
            structures_dir=None,  # Missing!
            bulk_name=None,       # Missing!
            metal_name=None,
            oxygen_name=None,
            miller_indices=None,
            input_slabs=None,
        )
    assert 'Surface thermodynamics' in str(exc_info.value)
```

## Backward Compatibility

The extraction is fully backward compatible:
- Same validation logic
- Same error messages
- Same logging output
- Returns the `needs_bulk` boolean for use by downstream code
