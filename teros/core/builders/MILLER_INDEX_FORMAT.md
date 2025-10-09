# Miller Index Format - Important Note

## Issue Fixed

The slab generation code expects a **single Miller index** as a list of three integers, not a list of tuples.

### ❌ INCORRECT Format (was causing error):
```python
miller_indices = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]  # List of tuples - ERROR!
```

**Error**: `TypeError: 'list' object cannot be interpreted as an integer`

### ✅ CORRECT Format (now used):
```python
miller_indices = [1, 0, 0]  # Single surface - (100)
```

## Current Limitation

**Only ONE surface can be generated at a time.**

To generate different surfaces, you need to run separate calculations:

```python
from teros.core.builders.default_ag3po4_builders import get_ag3po4_defaults
from teros.core.workgraph import build_core_workgraph_with_map

# Run 1: Generate (100) surface
defaults_100 = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    miller_indices=[1, 0, 0]  # (100) surface
)
wg_100 = build_core_workgraph_with_map(**defaults_100, name="Ag3PO4_100")

# Run 2: Generate (110) surface
defaults_110 = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    miller_indices=[1, 1, 0]  # (110) surface
)
wg_110 = build_core_workgraph_with_map(**defaults_110, name="Ag3PO4_110")

# Run 3: Generate (111) surface
defaults_111 = get_ag3po4_defaults(
    structures_dir="/path/to/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    miller_indices=[1, 1, 1]  # (111) surface
)
wg_111 = build_core_workgraph_with_map(**defaults_111, name="Ag3PO4_111")
```

## Default Miller Index

The default is the **(100) surface**:

```python
defaults = get_ag3po4_defaults(
    structures_dir="/path",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE"
)

print(defaults['miller_indices'])  # [1, 0, 0]
```

## Overriding Miller Index

To generate a different surface, override the `miller_indices` parameter:

### Option 1: At creation
```python
defaults = get_ag3po4_defaults(
    structures_dir="/path",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    miller_indices=[1, 1, 0]  # (110) surface
)
```

### Option 2: After creation
```python
defaults = get_ag3po4_defaults(
    structures_dir="/path",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE"
)
defaults['miller_indices'] = [1, 1, 1]  # Change to (111) surface
```

## Common Miller Indices for Ag3PO4

```python
# (100) surface - Default
miller_indices = [1, 0, 0]

# (110) surface
miller_indices = [1, 1, 0]

# (111) surface
miller_indices = [1, 1, 1]

# (210) surface
miller_indices = [2, 1, 0]
```

## Complete Example

```python
#!/usr/bin/env python
"""Generate and relax (100) surface of Ag3PO4."""

from aiida import load_profile
from teros.core.builders.default_ag3po4_builders import get_ag3po4_defaults
from teros.core.workgraph import build_core_workgraph_with_map

# Load AiiDA
load_profile(profile='psteros')

# Get defaults for (100) surface
defaults = get_ag3po4_defaults(
    structures_dir="/home/thiagotd/git/PS-TEROS/examples/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    miller_indices=[1, 0, 0]  # (100) surface
)

# Create and submit workgraph
wg = build_core_workgraph_with_map(
    **defaults,
    name="Ag3PO4_100_surface"
)
wg.submit(wait=False)

print(f"Submitted WorkGraph PK: {wg.pk}")
print(f"Surface: ({defaults['miller_indices'][0]} {defaults['miller_indices'][1]} {defaults['miller_indices'][2]})")
```

## Files Updated

All files have been updated with the correct format:

1. **`teros/core/builders/default_ag3po4_builders.py`**
   - Changed: `miller_indices = [1, 0, 0]`
   - Added documentation note about single surface limitation

2. **`examples/default_builders/test_default_builders.py`**
   - Updated assertion: `assert defaults['miller_indices'] == [1, 0, 0]`
   - All tests passing ✅

3. **`examples/default_builders/slabs_autogen_ag3po4.py`**
   - Updated examples to show correct format
   - Added note about single surface limitation

## Verification

Run tests to verify the format is correct:

```bash
cd /home/thiagotd/git/worktree/PS-TEROS/default-builders
source ~/envs/psteros/bin/activate
python examples/default_builders/test_default_builders.py
```

**Result**: ✅ All 5/5 tests passed

Check the Miller index value:

```bash
python -c "
from teros.core.builders.default_ag3po4_builders import get_ag3po4_defaults
d = get_ag3po4_defaults(structures_dir='/tmp', code_label='test', potential_family='PBE')
print(f'Miller index: {d[\"miller_indices\"]}')
print(f'Type: {type(d[\"miller_indices\"])}')
"
```

**Output**:
```
Miller index: [1, 0, 0]
Type: <class 'list'>
```

## Summary

- ✅ **Fixed**: Miller index format now correct: `[h, k, l]`
- ✅ **Default**: `[1, 0, 0]` - (100) surface
- ✅ **Tests**: All passing (5/5)
- ✅ **Documentation**: Updated with correct examples
- ⚠️ **Limitation**: Only one surface per calculation
- ✅ **Workaround**: Run multiple calculations for multiple surfaces

---

**Updated**: 2025-01-09
**Status**: ✅ Fixed and tested
**Module**: `teros.core.builders.default_ag3po4_builders`
