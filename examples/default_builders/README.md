# Default Builders Examples

This directory contains examples demonstrating the use of the `teros.default_builders` module.

## Files

### `slabs_input_relax_ag3po4.py` (Original)
The original example file with all VASP parameters explicitly defined (~380 lines).

**Use this to**: Understand what all the default parameters are doing.

### `slabs_input_relax_ag3po4_simplified.py` (Simplified)
Simplified example using the `default_builders` module (~180 lines).

**Use this to**: See how to use default builders in your own workflows.

**Key differences:**
- ~200 lines of parameter definitions replaced with one function call
- Cleaner, more maintainable code
- Easy to override specific parameters
- Same functionality as the original

### `test_default_builders.py`
Test script that validates the default_builders module functionality.

**Use this to**: 
- Verify the module works correctly after modifications
- Understand how the module handles parameter overrides
- See examples of all module features

## Quick Start

### 1. Test the Module

```bash
cd /home/thiagotd/git/worktree/PS-TEROS/default-builders
source ~/envs/psteros/bin/activate
python examples/default_builders/test_default_builders.py
```

### 2. Run the Simplified Example

```bash
# Make sure your slab files exist
ls /home/thiagotd/git/PS-TEROS/examples/slabs/input_structures/

# Run the workflow
source ~/envs/psteros/bin/activate
python examples/default_builders/slabs_input_relax_ag3po4_simplified.py
```

### 3. Create Your Own Workflow

Copy the simplified example and modify:

```python
from teros.default_builders import get_ag3po4_defaults
from teros.core.workgraph import build_core_workgraph_with_map

# Get defaults
defaults = get_ag3po4_defaults(
    structures_dir="/your/path/to/structures",
    code_label="your-vasp-code",
    potential_family="PBE",
    # Add your overrides here
    bulk_parameters={'ENCUT': 600}
)

# Create workgraph
wg = build_core_workgraph_with_map(
    **defaults,
    input_slabs=your_slabs,
    name="YourCalculation"
)
```

## Comparison

### Before (Original)

```python
# ~200 lines of this:
bulk_parameters = {
    "PREC": "Accurate",
    "ENCUT": 520,
    "EDIFF": 1e-6,
    # ... 10 more parameters
}

bulk_options = {
    "resources": {
        "num_machines": 1,
        # ... more options
    },
    # ... more options
}

metal_parameters = {
    # ... another 15 parameters
}

# ... and so on for metal, nonmetal, oxygen, slab
```

### After (Simplified)

```python
# Just 1 line:
defaults = get_ag3po4_defaults(
    structures_dir="/path",
    code_label="vasp-code",
    potential_family="PBE"
)

# Optional overrides:
defaults['bulk_parameters']['ENCUT'] = 600
```

## Documentation

See the full documentation: `/home/thiagotd/git/worktree/PS-TEROS/default-builders/docs/default_builders_guide.md`

## Benefits

1. **Less code**: ~50% reduction in lines of code
2. **Fewer errors**: Default values are tested and validated
3. **Consistency**: Same defaults across all calculations
4. **Maintainability**: Update defaults in one place
5. **Flexibility**: Easy to override any parameter
6. **Readability**: Intent is clear, not buried in parameters
