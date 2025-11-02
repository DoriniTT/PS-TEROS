# Custom VASP Calculation Module Design

**Date:** 2025-11-01
**Branch:** `feature-custom-calculation`
**Status:** Design approved, ready for implementation

## Purpose

This module runs arbitrary VASP calculations through PS-TEROS WorkGraphs. Users define all VASP parameters in their run scripts and submit calculations on single structures or multiple structures with full control over builder inputs.

## Requirements

1. **Flexibility**: Support any VASP calculation type (relaxation, static, DOS, custom INCAR settings)
2. **Control**: User defines all builder inputs in the run script
3. **Scope**: Handle single structure or multiple structures
4. **Output**: Extract total energy, relaxed structure, and raw VASP outputs
5. **Integration**: Follow existing PS-TEROS patterns (like `surface_hydroxylation` module)

## Module Structure

```
teros/core/custom_calculation/
├── __init__.py          # Exports: build_custom_calculation_workgraph, get_custom_results
├── workgraph.py         # Main builder function and result extraction
└── tasks.py             # Helper tasks: extract_energy, extract_structure
```

The module contains three files:
- `workgraph.py`: Creates WorkGraph with VaspTask(s)
- `tasks.py`: Extracts energy, structure, and files from VASP output
- `__init__.py`: Exports public API

## API Design

### Main Function

```python
build_custom_calculation_workgraph(
    structure,          # StructureData or list[StructureData]
    code_label,         # str: VASP code label
    builder_inputs,     # dict or list[dict]: full VASP builder parameters
    name='custom_calc'  # str: WorkGraph name
) -> WorkGraph
```

### Single Structure Example

```python
from teros.core.custom_calculation import build_custom_calculation_workgraph

# Define all builder inputs
builder_inputs = {
    'parameters': {'incar': {
        'PREC': 'Accurate',
        'ENCUT': 500,
        'EDIFF': 1e-5,
        'ISMEAR': 0,
        'SIGMA': 0.05,
    }},
    'options': {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
    },
    'kpoints_spacing': 0.3,
    'potential_family': 'PBE.54',
    'potential_mapping': {'Ag': 'Ag', 'O': 'O'},
    'clean_workdir': True,
}

# Build and submit
wg = build_custom_calculation_workgraph(
    structure=my_structure,
    code_label='VASP-6.4.1@cluster02',
    builder_inputs=builder_inputs,
    name='my_custom_calc'
)

wg.submit(wait=True)
```

### Multiple Structures - Same Settings

```python
# Same builder inputs for all structures
wg = build_custom_calculation_workgraph(
    structure=[struct1, struct2, struct3],
    code_label='VASP-6.4.1@cluster02',
    builder_inputs=builder_inputs,  # Single dict reused
)
```

### Multiple Structures - Different Settings

```python
# Different builder inputs for each structure
builder_inputs_list = [
    {...},  # settings for structure 1
    {...},  # settings for structure 2
    {...},  # settings for structure 3
]

wg = build_custom_calculation_workgraph(
    structure=[struct1, struct2, struct3],
    code_label='VASP-6.4.1@cluster02',
    builder_inputs=builder_inputs_list,  # List matches structure list
)
```

## WorkGraph Structure

### Single Structure Flow

```
Input: structure, code, builder_inputs
  ↓
VaspTask (runs VASP calculation)
  ↓
ExtractEnergyTask → energy (Float)
ExtractStructureTask → structure (StructureData)
  ↓
WorkGraph outputs:
  - energy: Float
  - structure: StructureData
  - misc: Dict (raw VASP outputs)
```

### Multiple Structures Flow

```
Input: [struct1, struct2, struct3], code, [builder1, builder2, builder3]
  ↓
Map over structures → VaspTask for each
  ↓
Extract tasks for each → parallel extraction
  ↓
WorkGraph outputs:
  - energies: [Float, Float, Float]
  - structures: [StructureData, StructureData, StructureData]
  - misc: [Dict, Dict, Dict]
```

## Result Extraction

### Direct Access

```python
wg.submit(wait=True)

# Single structure
energy = wg.outputs.energy.value
structure = wg.outputs.structure
misc = wg.outputs.misc.get_dict()

# Multiple structures
energies = [e.value for e in wg.outputs.energies]
structures = wg.outputs.structures
misc_list = [m.get_dict() for m in wg.outputs.misc]
```

### Helper Function

```python
from teros.core.custom_calculation import get_custom_results

results = get_custom_results(wg)
# Returns:
# {
#     'energies': [value1, value2, ...] or single value,
#     'structures': [struct1, struct2, ...] or single structure,
#     'misc': [dict1, dict2, ...] or single dict
# }
```

## Error Handling

- Failed VaspTask shows `EXCEPTED` or `FAILED` state in WorkGraph
- Use standard AiiDA debugging: `verdi process show <PK>`, `verdi process report <PK>`
- Module adds no extra error handling—keeps implementation simple and transparent
- User checks WorkGraph state before accessing results

## Implementation Details

### Builder Input Structure

The `builder_inputs` dict contains standard VASP parameters:

```python
{
    'parameters': {'incar': {...}},     # INCAR settings
    'options': {...},                    # Compute resources
    'kpoints_spacing': float,            # K-points spacing
    'potential_family': str,             # Pseudopotential family
    'potential_mapping': dict,           # Element to potential mapping
    'clean_workdir': bool,               # Clean work directory after run
    'settings': orm.Dict (optional),     # Additional VASP settings
}
```

### Auto-Detection Logic

The builder function detects input type:
- If `structure` is StructureData → create one VaspTask
- If `structure` is list → create VaspTask for each structure
- If `builder_inputs` is dict → reuse for all structures
- If `builder_inputs` is list → must match structure list length

### Task Implementation

**VaspTask creation:**
```python
from aiida.plugins import WorkflowFactory
from aiida_workgraph import task

VaspWorkChain = WorkflowFactory('vasp.v2.vasp')
VaspTask = task(VaspWorkChain)
```

**Energy extraction:**
```python
def extract_energy(misc):
    """Extract total energy from VASP misc output."""
    return misc['total_energies']['energy_extrapolated']
```

**Structure extraction:**
```python
def extract_structure(misc):
    """Extract relaxed structure from VASP misc output."""
    return misc['structure']
```

## Testing Strategy

Create test examples in `examples/custom_calculation/`:

1. **Single structure test** (`test_single.py`):
   - Load one structure
   - Define builder inputs
   - Run calculation
   - Verify outputs (energy, structure, misc)
   - Check WorkGraph completes with `[0]` status

2. **Multiple structures test** (`test_multiple.py`):
   - Load three structures
   - Use same builder inputs for all
   - Run calculations in parallel
   - Verify all outputs
   - Check WorkGraph completes successfully

3. **Different settings test** (`test_different_settings.py`):
   - Load three structures
   - Define different builder inputs for each
   - Run calculations
   - Verify different INCAR settings applied correctly

Test scripts follow `examples/hydroxylation_with_bulk_reference/` pattern.

## Integration with PS-TEROS

### Follows Existing Patterns

1. **Module structure**: Same as `surface_hydroxylation` (separate folder with `__init__.py`, `workgraph.py`, `tasks.py`)
2. **Import style**: `from teros.core.custom_calculation import build_custom_calculation_workgraph`
3. **Builder parameters**: Match existing VASP patterns in `workgraph.py` and `hydroxylation`
4. **Testing location**: `examples/custom_calculation/` with complete test cases
5. **Code loading**: User specifies code label in run script

### No Breaking Changes

This module:
- Adds new functionality only
- Does not modify existing modules
- Can be imported and used alongside other modules
- Requires no changes to existing code

## File Checklist

Implementation requires these files:

```
teros/core/custom_calculation/
├── __init__.py
├── workgraph.py
└── tasks.py

examples/custom_calculation/
├── test_single.py
├── test_multiple.py
└── test_different_settings.py
```

## Success Criteria

Implementation succeeds when:

1. All test scripts run without errors
2. Main WorkGraph node returns `[0]` (successful completion)
3. Outputs appear correctly: `energy`, `structure`, `misc` for single structure
4. Outputs appear correctly: `energies`, `structures`, `misc` for multiple structures
5. User has full control over all VASP parameters in run script
6. Module follows PS-TEROS conventions and patterns

## Next Steps

1. Set up git worktree for `feature-custom-calculation` branch
2. Create implementation plan with detailed tasks
3. Implement module following TDD (test-driven development)
4. Test with production-like parameters
5. Merge to `develop` after successful validation
