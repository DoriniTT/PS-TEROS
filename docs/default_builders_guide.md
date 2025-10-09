# Default Builders Guide

## Overview

The `teros.default_builders` module provides pre-configured VASP calculation parameters, scheduler options, and structure configurations for specific material systems. This reduces boilerplate code and ensures consistency across calculations while still allowing easy customization.

## Features

- **Pre-configured defaults**: All VASP INCAR parameters, scheduler options, and structure files for supported materials
- **Easy customization**: Override any parameter with simple dictionary syntax
- **Deep merging**: Partial overrides merge with defaults instead of replacing them
- **Type safety**: Uses the same structure as `build_core_workgraph_with_map()`
- **Extensible**: Easy to add new materials

## Supported Materials

### Ag₃PO₄ (Silver Phosphate)

Function: `get_ag3po4_defaults()`

Includes defaults for:
- Bulk Ag₃PO₄ relaxation
- Metal (Ag) reference
- Nonmetal (P) reference
- Oxygen (O₂) reference
- Slab relaxation with dipole corrections

## Basic Usage

### Simple Example (No Overrides)

```python
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph_with_map
from teros.default_builders import get_ag3po4_defaults

# Load AiiDA profile
load_profile(profile='psteros')

# Get all default parameters for Ag3PO4
defaults = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE"
)

# Create workgraph with defaults
wg = build_core_workgraph_with_map(
    **defaults,  # Unpack all defaults
    input_slabs=my_input_slabs,  # Add run-specific parameters
    name="MyCalculation"
)
```

### Advanced Example (With Overrides)

```python
from teros.default_builders import get_ag3po4_defaults

# Get defaults with custom overrides
defaults = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    # Override specific VASP parameters
    bulk_parameters={'ENCUT': 600, 'EDIFF': 1e-7},
    slab_parameters={'EDIFFG': -0.01},
    # Override scheduler options
    slab_options={
        'resources': {
            'num_machines': 2,
            'num_cores_per_machine': 40
        }
    },
    # Override k-points
    kpoints_spacing=0.25,
    slab_kpoints_spacing=0.25
)

# Create workgraph
wg = build_core_workgraph_with_map(**defaults, input_slabs=my_slabs)
```

### Manual Override After Getting Defaults

```python
from teros.default_builders import get_ag3po4_defaults, update_builder_params

# Get base defaults
defaults = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE"
)

# Apply overrides using helper function
defaults = update_builder_params(defaults, {
    'bulk_parameters': {'ENCUT': 600},
    'slab_options': {
        'resources': {'num_machines': 2}
    }
})

# Or modify directly (shallow override)
defaults['bulk_parameters']['ENCUT'] = 600
defaults['thermodynamics_sampling'] = 200

# Create workgraph
wg = build_core_workgraph_with_map(**defaults, input_slabs=my_slabs)
```

## API Reference

### `get_ag3po4_defaults()`

Get default builder parameters for Ag₃PO₄ material system.

**Arguments:**
- `structures_dir` (str, optional): Path to directory containing structure files
- `code_label` (str, optional): VASP code label (e.g., "VASP-VTST-6.4.3@bohr")
- `potential_family` (str, optional): Potential family name (e.g., "PBE")
- `**overrides`: Any parameter to override (deep merged with defaults)

**Returns:**
- `dict`: Complete builder parameters dictionary compatible with `build_core_workgraph_with_map()`

**Default Parameters:**

Bulk relaxation:
```python
{
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 100,
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}
```

Metal (Ag) reference:
```python
{
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 100,
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}
```

Nonmetal (P) reference:
```python
{
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 100,
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
}
```

Oxygen (O₂) reference:
```python
{
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.01,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 100,
    "EDIFFG": -0.1,
    "ALGO": "Normal",
    "LREAL": False,
    "LWAVE": False,
    "LCHARG": False,
}
```

Slab relaxation:
```python
{
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 100,
    "EDIFFG": -0.02,
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
    "IDIPOL": 3,
    "LDIPOL": True,
    "LVHAR": True,
}
```

### `update_builder_params()`

Deep merge override parameters into default parameters.

**Arguments:**
- `defaults` (dict): Default parameters dictionary
- `overrides` (dict): Override parameters to merge into defaults

**Returns:**
- `dict`: Merged parameters dictionary

**Example:**
```python
defaults = {'params': {'ENCUT': 520, 'EDIFF': 1e-6}}
overrides = {'params': {'ENCUT': 600}}
result = update_builder_params(defaults, overrides)
# result = {'params': {'ENCUT': 600, 'EDIFF': 1e-6}}
```

## Complete Example

See the simplified example file:
```
examples/default_builders/slabs_input_relax_ag3po4_simplified.py
```

This file demonstrates how to use the default builders module and compares to the original full example (which had ~200 lines of parameter definitions replaced with a single function call).

## Adding New Materials

To add defaults for a new material, follow this template:

```python
def get_<material>_defaults(
    structures_dir=None,
    code_label=None,
    potential_family=None,
    **overrides
):
    """
    Get default builder parameters for <Material> material system.
    
    Args:
        structures_dir (str, optional): Path to directory containing structure files.
        code_label (str, optional): VASP code label.
        potential_family (str, optional): Potential family name.
        **overrides: Any parameter to override.
        
    Returns:
        dict: Complete builder parameters dictionary.
    """
    # Define bulk_parameters
    bulk_parameters = {...}
    
    # Define bulk_options
    bulk_options = {...}
    
    # Define metal_parameters, metal_options
    # Define nonmetal_parameters, nonmetal_options
    # Define oxygen_parameters, oxygen_options
    # Define slab_parameters, slab_options
    
    # Define potential mappings
    bulk_potential_mapping = {...}
    # ... etc
    
    # Define structure file names
    bulk_name = "material.cif"
    metal_name = "metal.cif"
    # ... etc
    
    # Build defaults dictionary
    defaults = {
        "structures_dir": structures_dir,
        "bulk_name": bulk_name,
        # ... all other parameters
    }
    
    # Apply overrides
    if overrides:
        defaults = update_builder_params(defaults, overrides)
    
    return defaults
```

## Best Practices

1. **Start with defaults**: Always start with the material defaults and only override what you need
2. **Use deep merging**: Use `update_builder_params()` for complex overrides to preserve nested structure
3. **Document overrides**: Add comments explaining why you're overriding specific parameters
4. **Test incrementally**: Test with defaults first, then add overrides one at a time
5. **Check consistency**: Ensure overridden parameters are physically meaningful for your system

## Comparison: Before and After

### Before (Original File)

~380 lines with:
- ~200 lines of parameter definitions
- Repetitive code for each calculation type
- Hard to maintain and modify
- Error-prone when copying parameters

### After (With Default Builders)

~180 lines with:
- ~3 lines to get defaults
- Optional parameter overrides only where needed
- Clean and maintainable
- Easy to understand intent

**Code reduction: ~50%** while maintaining full flexibility!

## Troubleshooting

### Issue: Parameters not being overridden

**Problem**: Using assignment (`=`) on nested dictionaries replaces the entire dict instead of merging.

**Solution**: Use `update_builder_params()` for deep merging:
```python
# Wrong: Replaces entire dictionary
defaults['bulk_parameters'] = {'ENCUT': 600}  # Loses all other parameters!

# Right: Deep merges parameters
defaults = update_builder_params(defaults, {
    'bulk_parameters': {'ENCUT': 600}  # Keeps other parameters
})
```

### Issue: Missing required parameters

**Problem**: `structures_dir`, `code_label`, or `potential_family` is None.

**Solution**: Always provide these required parameters:
```python
defaults = get_ag3po4_defaults(
    structures_dir="/path/to/structures",  # Required
    code_label="VASP-VTST-6.4.3@bohr",     # Required
    potential_family="PBE"                  # Required
)
```

### Issue: Structure files not found

**Problem**: Structure files (ag3po4.cif, Ag.cif, P.cif, O2.cif) don't exist in `structures_dir`.

**Solution**: Ensure all required structure files exist:
```bash
ls /path/to/structures/
# Should show: ag3po4.cif, Ag.cif, P.cif, O2.cif
```

Or override the structure names:
```python
defaults = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    bulk_name="my_ag3po4.vasp",
    metal_name="my_Ag.cif",
    # ... etc
)
```
