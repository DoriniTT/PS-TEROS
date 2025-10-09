# User-Provided Slabs Feature

## Summary

PS-TEROS now supports using pre-generated slab structures as input, giving users full flexibility to provide their own slabs instead of relying on automatic slab generation. This feature maintains backward compatibility while adding new capabilities.

## Feature Overview

### Previous Behavior
- Slabs were always generated automatically from the relaxed bulk structure
- Required: `miller_indices`, `min_slab_thickness`, `min_vacuum_thickness`
- Used Pymatgen's SlabGenerator internally

### New Behavior
- **Option 1 (default)**: Continue using automatic generation (backward compatible)
- **Option 2 (new)**: Provide pre-generated slab structures via `input_slabs` parameter
- When `input_slabs` is provided, slab generation is skipped entirely

## API Changes

### Modified Functions

All functions in `teros/core/workgraph.py` now accept an optional `input_slabs` parameter:

1. `core_workgraph()` - Main task graph
2. `build_core_workgraph()` - Builder function
3. `build_core_workgraph_with_map()` - Alias for scatter-gather pattern

### New Parameter

```python
input_slabs: dict = None
```

- **Type**: Dictionary mapping string keys to `orm.StructureData` nodes
- **Format**: `{"term_0": StructureData, "term_1": StructureData, ...}`
- **Default**: `None` (automatic generation mode)

### Parameter Changes

When `input_slabs` is provided:
- `miller_indices` - Optional (ignored)
- `min_slab_thickness` - Optional (ignored)
- `min_vacuum_thickness` - Optional (ignored)
- Other generation parameters (lll_reduce, etc.) - Optional (ignored)

When `input_slabs` is `None`:
- `miller_indices` - **Required**
- `min_slab_thickness` - **Required**
- `min_vacuum_thickness` - **Required**

## Implementation Details

### Code Location
- **File**: `teros/core/workgraph.py`
- **Lines**: ~244-267 (slab generation section)

### Logic Flow
```python
if input_slabs is not None:
    # Use user-provided slabs
    slab_namespace = input_slabs
else:
    # Generate slabs automatically (original behavior)
    if miller_indices is None or min_slab_thickness is None or min_vacuum_thickness is None:
        raise ValueError("Generation parameters required when input_slabs not provided")
    slab_namespace = generate_slab_structures(...)
```

## Usage Examples

### Example 1: Automatic Generation (Original)

```python
wg = build_core_workgraph_with_map(
    structures_dir=structures_dir,
    bulk_name="ag3po4.cif",
    # ... other params ...
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    relax_slabs=True,
)
```

### Example 2: User-Provided Slabs (New)

```python
from ase.io import read
from aiida import orm

# Load slabs
input_slabs = {}
for idx, file in enumerate(["slab_0.cif", "slab_1.cif"]):
    atoms = read(file)
    input_slabs[f"term_{idx}"] = orm.StructureData(ase=atoms)

# Build workflow
wg = build_core_workgraph_with_map(
    structures_dir=structures_dir,
    bulk_name="ag3po4.cif",
    # ... other params ...
    input_slabs=input_slabs,  # Provide slabs
    relax_slabs=True,
)
```

## Files Modified

### Core Files
1. **`teros/core/workgraph.py`**
   - Added `input_slabs` parameter to `core_workgraph()`
   - Added conditional logic for slab generation vs. user input
   - Updated docstrings
   - Made generation parameters optional
   - Added `input_slabs` to `build_core_workgraph()`
   - Added `input_slabs` to `build_core_workgraph_with_map()`

### Example Files Created
1. **`examples/slabs/slabs_input_relax.py`**
   - Complete working example using user-provided slabs
   - Based on `examples/backup/slabs/slabs_relax.py`
   - Shows how to load structures from files
   - Demonstrates workflow setup and submission

2. **`examples/slabs/compare_modes.py`**
   - Comparison of automatic vs. user-provided modes
   - Demonstrates API differences
   - Shows use cases for each mode

### Documentation Created
1. **`examples/slabs/README_INPUT_SLABS.md`**
   - Comprehensive guide for using user-provided slabs
   - Usage instructions and best practices
   - Troubleshooting tips

2. **`examples/slabs/input_structures/README.md`**
   - Instructions for the input structures directory
   - File format and naming conventions

3. **`docs/USER_PROVIDED_SLABS.md`** (this file)
   - Technical documentation
   - API changes and implementation details

## Backward Compatibility

✅ **Fully backward compatible**
- Existing scripts continue to work without modification
- Default behavior unchanged (automatic generation)
- All original parameters function as before
- No breaking changes to API

## Testing

### Syntax Check
```bash
python -m py_compile teros/core/workgraph.py
python -m py_compile examples/slabs/slabs_input_relax.py
```

### Validation
- Comparison script runs successfully
- No syntax errors in modified files
- Import statements work correctly

## Benefits

1. **Flexibility**: Use any slab generation method
2. **Control**: Exact control over surface structures
3. **Reproducibility**: Use structures from literature
4. **Efficiency**: Skip generation when structures are available
5. **Compatibility**: Works with all ASE-supported formats
6. **Integration**: Seamless integration with existing workflows

## Use Cases

- Reproducing published structures
- Custom surface reconstructions
- Adding adsorbates or defects
- Non-standard surface orientations
- Manual slab editing
- Using specialized generation tools
- Testing specific configurations

## Future Enhancements

Possible future improvements:
- Mixed mode (some generated, some provided)
- Automatic validation of slab structures
- Helper functions for common slab operations
- Integration with other slab generation tools

## Migration Guide

No migration needed! Existing code works as-is. To use new feature:

1. Load your slab structures:
   ```python
   input_slabs = {"term_0": structure1, "term_1": structure2}
   ```

2. Add to workflow builder:
   ```python
   wg = build_core_workgraph_with_map(..., input_slabs=input_slabs)
   ```

3. Remove generation parameters (optional):
   - Can remove `miller_indices`, `min_slab_thickness`, etc.
   - Or leave them in (will be ignored)

## Support

For questions or issues:
- See examples in `examples/slabs/`
- Read guide in `examples/slabs/README_INPUT_SLABS.md`
- Check comparison in `examples/slabs/compare_modes.py`
