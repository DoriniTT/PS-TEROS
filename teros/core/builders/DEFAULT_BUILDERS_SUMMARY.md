# Default Builders Module - Summary

## Overview

A new module `teros.default_builders` has been created to store default VASP calculation parameters for specific material systems, reducing boilerplate code and improving maintainability.

## What Was Created

### 1. Core Module
**File**: `teros/default_builders.py`

Contains:
- `get_ag3po4_defaults()`: Returns all default parameters for Ag₃PO₄ system
- `update_builder_params()`: Helper for deep merging parameter overrides
- Template for adding new materials

### 2. Documentation
**File**: `docs/default_builders_guide.md`

Comprehensive guide including:
- API reference
- Usage examples
- Default parameter values
- Troubleshooting guide
- Best practices

### 3. Examples
**Directory**: `examples/default_builders/`

**Files**:
- `slabs_input_relax_ag3po4.py` - Original example (reference)
- `slabs_input_relax_ag3po4_simplified.py` - Simplified example using default_builders
- `test_default_builders.py` - Test suite for the module
- `README.md` - Quick start guide

## Key Features

### ✨ Simplification
Reduces ~200 lines of parameter definitions to a single function call:

```python
# Before: 200 lines of parameter definitions
bulk_parameters = {...}
bulk_options = {...}
metal_parameters = {...}
# ... etc

# After: 1 function call
defaults = get_ag3po4_defaults(
    structures_dir="/path",
    code_label="vasp-code",
    potential_family="PBE"
)
```

### 🔧 Easy Customization
Three ways to override parameters:

1. **At creation**:
```python
defaults = get_ag3po4_defaults(
    structures_dir="/path",
    code_label="vasp-code",
    potential_family="PBE",
    bulk_parameters={'ENCUT': 600}  # Override
)
```

2. **Using helper function**:
```python
defaults = update_builder_params(defaults, {
    'bulk_parameters': {'ENCUT': 600}
})
```

3. **Direct modification**:
```python
defaults['bulk_parameters']['ENCUT'] = 600
```

### 🏗️ Deep Merging
Overrides merge with defaults instead of replacing them:
```python
# Override ENCUT but keep all other parameters
defaults = get_ag3po4_defaults(
    bulk_parameters={'ENCUT': 600}
)
# Result: ENCUT=600, but EDIFF, ISMEAR, etc. are preserved
```

### 📦 Complete Integration
Works seamlessly with existing PS-TEROS code:
```python
wg = build_core_workgraph_with_map(
    **defaults,  # Unpack all defaults
    input_slabs=my_slabs,  # Add run-specific parameters
    name="MyCalculation"
)
```

## Default Parameters Included

For Ag₃PO₄ system:

### Bulk Relaxation
- VASP INCAR parameters (ENCUT=520, EDIFF=1e-6, ISMEAR=0, etc.)
- Scheduler options (40 cores, par40 queue)
- Structure file: `ag3po4.cif`
- Potential mapping: Ag→Ag, P→P, O→O

### Metal Reference (Ag)
- VASP INCAR parameters (ISMEAR=1, SIGMA=0.2, etc.)
- Scheduler options
- Structure file: `Ag.cif`
- Potential mapping: Ag→Ag

### Nonmetal Reference (P)
- VASP INCAR parameters (ISMEAR=0, SIGMA=0.05, etc.)
- Scheduler options
- Structure file: `P.cif`
- Potential mapping: P→P

### Oxygen Reference (O₂)
- VASP INCAR parameters (ISMEAR=0, SIGMA=0.01, ISIF=2, etc.)
- Scheduler options
- Structure file: `O2.cif`
- Potential mapping: O→O

### Slab Relaxation
- VASP INCAR parameters (EDIFFG=-0.02, IDIPOL=3, LDIPOL=True, etc.)
- Scheduler options
- Potential mapping: Ag→Ag, P→P, O→O
- K-points spacing: 0.3

### Other Settings
- K-points spacing: 0.3
- Clean workdir: True
- Relax slabs: True
- Compute thermodynamics: True
- Thermodynamics sampling: 100

## Testing

All functionality is tested and validated:

```bash
cd /home/thiagotd/git/worktree/PS-TEROS/default-builders
source ~/envs/psteros/bin/activate
python examples/default_builders/test_default_builders.py
```

**Result**: ✓ All 5 tests passed

Tests validate:
1. Default parameters load correctly
2. Parameters can be overridden at creation
3. Deep merging works correctly
4. Structure information is correct
5. Dictionary unpacking works for function calls

## Usage Example

### Complete Workflow

```python
#!/home/thiagotd/envs/psteros/bin/python
from aiida import load_profile, orm
from ase.io import read
from teros.core.workgraph import build_core_workgraph_with_map
from teros.default_builders import get_ag3po4_defaults

# Load AiiDA
load_profile(profile='psteros')

# Load your slabs
input_slabs = {}
for i in range(3):
    atoms = read(f"slab_term_{i}.cif")
    input_slabs[f"term_{i}"] = orm.StructureData(ase=atoms)

# Get defaults with optional overrides
defaults = get_ag3po4_defaults(
    structures_dir="/home/thiagotd/git/PS-TEROS/examples/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    # Optional: override specific parameters
    bulk_parameters={'ENCUT': 600},
    slab_parameters={'EDIFFG': -0.01}
)

# Create and submit workgraph
wg = build_core_workgraph_with_map(
    **defaults,
    input_slabs=input_slabs,
    name="Ag3PO4_Slabs"
)

wg.submit(wait=False)
print(f"Submitted WorkGraph PK: {wg.pk}")
```

## Benefits

### For Users
- **Less typing**: ~50% code reduction
- **Fewer errors**: Validated defaults
- **Faster development**: Focus on science, not parameters
- **Easier learning**: Clear examples

### For Maintainers
- **Single source of truth**: Update defaults in one place
- **Version control**: Track parameter changes
- **Testing**: Validate defaults automatically
- **Documentation**: Parameters are documented

### For Collaboration
- **Consistency**: Everyone uses same defaults
- **Reproducibility**: Easy to share parameters
- **Transparency**: Clear what's being used

## Next Steps

### Using the Module
1. Review the documentation: `docs/default_builders_guide.md`
2. Run the test: `examples/default_builders/test_default_builders.py`
3. Try the simplified example: `examples/default_builders/slabs_input_relax_ag3po4_simplified.py`
4. Create your own workflow using `get_ag3po4_defaults()`

### Extending the Module
To add new materials, follow the template in `teros/default_builders.py`:

```python
def get_<material>_defaults(
    structures_dir=None,
    code_label=None,
    potential_family=None,
    **overrides
):
    """Get default builder parameters for <Material>."""
    # Define all parameters
    # Build defaults dictionary
    # Apply overrides
    return defaults
```

## Files Modified/Created

```
teros/
  └── default_builders.py (NEW)         - Core module

docs/
  └── default_builders_guide.md (NEW)  - Full documentation

examples/default_builders/
  ├── README.md (NEW)                   - Quick start
  ├── slabs_input_relax_ag3po4.py      - Original (unchanged)
  ├── slabs_input_relax_ag3po4_simplified.py (NEW) - Using module
  └── test_default_builders.py (NEW)   - Tests

DEFAULT_BUILDERS_SUMMARY.md (NEW)      - This file
```

## Important Notes

1. **After modifications**: Run `verdi daemon restart` to reload changes
2. **Structure files**: Ensure all structure files exist in `structures_dir`
3. **Overrides**: Use `update_builder_params()` for deep merging
4. **Testing**: Run tests after any module modifications

## Support

- **Documentation**: `docs/default_builders_guide.md`
- **Examples**: `examples/default_builders/`
- **Tests**: `examples/default_builders/test_default_builders.py`

## Version

- **Created**: 2025
- **Module**: `teros.default_builders`
- **First material**: Ag₃PO₄
- **Status**: ✅ Tested and validated
