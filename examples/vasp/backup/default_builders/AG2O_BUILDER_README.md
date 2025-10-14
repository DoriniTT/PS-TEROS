# Ag2O Default Builder Documentation

## Overview

The `default_ag2o_builders` module provides pre-configured parameters for **binary oxide** Ag₂O calculations. This is specifically designed for binary metal oxides where only metal and oxygen references are physically meaningful.

## Key Differences from Ag3PO4

### Binary Oxide vs Ternary Oxide

**Ag₂O (Binary Oxide)**:
- Contains: Ag (metal) + O (oxygen) only
- References needed: Ag and O₂
- Nonmetal reference: **Dummy** (uses Ag parameters)
- Surface energy: γ(Δμ_O) - function of oxygen chemical potential
- Thermodynamics: Binary oxide model

**Ag₃PO₄ (Ternary Oxide)**:
- Contains: Ag (metal) + P (nonmetal) + O (oxygen)
- References needed: Ag, P, and O₂
- Nonmetal reference: **Physical** (P parameters)
- Surface energy: γ(Δμ_P, Δμ_O) - function of both chemical potentials
- Thermodynamics: Ternary oxide model

## Module Location

```
teros/core/builders/default_ag2o_builders.py
```

## Parameters Included

Total: **37 parameters** (6 more than Ag3PO4 due to additional slab generation options)

### Structure Files (4)
- `bulk_name`: "ag2o.cif"
- `metal_name`: "Ag.cif"
- `nonmetal_name`: "Ag.cif" ⚠️ **Dummy for binary oxide**
- `oxygen_name`: "O2.cif"

### Slab Generation (9)
- `miller_indices`: [1, 0, 0] (100) surface
- `min_slab_thickness`: 10.0 Å
- `min_vacuum_thickness`: 15.0 Å
- `lll_reduce`: True
- `center_slab`: True
- `symmetrize`: True
- `primitive`: True
- `in_unit_planes`: False
- `max_normal_search`: None

### Potential Mappings (5)
- `bulk_potential_mapping`: {"Ag": "Ag", "O": "O"}
- `metal_potential_mapping`: {"Ag": "Ag"}
- `nonmetal_potential_mapping`: {"Ag": "Ag"} ⚠️ **Dummy**
- `oxygen_potential_mapping`: {"O": "O"}
- `slab_potential_mapping`: {"Ag": "Ag", "O": "O"}

### Plus: All standard parameters (19)
- K-points settings (2)
- Bulk relaxation (2)
- Metal reference (2)
- Nonmetal reference (2) - dummy
- Oxygen reference (2)
- Slab relaxation (2)
- Other settings (7)

## Basic Usage

### Simple Example

```python
from aiida import load_profile
from teros.core.builders.default_ag2o_builders import get_ag2o_defaults
from teros.core.workgraph import build_core_workgraph_with_map

# Load AiiDA
load_profile(profile='psteros')

# Get defaults
defaults = get_ag2o_defaults(
    structures_dir="/home/thiagotd/git/PS-TEROS/examples/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE"
)

# Create workgraph - slab will be generated automatically
wg = build_core_workgraph_with_map(
    **defaults,
    name="Ag2O_100_surface"
)

wg.submit(wait=False)
print(f"Submitted WorkGraph PK: {wg.pk}")
```

### With Overrides

```python
# Get defaults with overrides
defaults = get_ag2o_defaults(
    structures_dir="/path/to/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    # Override Miller index
    miller_indices=[1, 1, 0],  # (110) surface
    # Override thickness
    min_slab_thickness=15.0,
    # Override symmetrization
    symmetrize=False,
    # Override VASP parameters
    bulk_parameters={'ENCUT': 600},
    slab_parameters={'EDIFFG': -0.01}
)

wg = build_core_workgraph_with_map(**defaults, name="Ag2O_110")
```

## Common Miller Indices for Ag2O

```python
# (100) surface - Default
miller_indices = [1, 0, 0]

# (110) surface
miller_indices = [1, 1, 0]

# (111) surface
miller_indices = [1, 1, 1]
```

## Important Notes

### 1. Binary Oxide System

Ag₂O is a **binary oxide** containing only metal (Ag) and oxygen (O). The nonmetal reference is **not physically meaningful** but is required by the PS-TEROS framework for compatibility.

**What this means:**
- Nonmetal parameters use the same values as metal (Ag)
- Nonmetal reference calculation will run but doesn't affect binary thermodynamics
- Surface energies are γ(Δμ_O) - only function of oxygen chemical potential

### 2. Miller Index Format

Only **ONE** Miller index at a time: `[h, k, l]` not `[(h, k, l), ...]`

```python
# ✅ CORRECT
miller_indices = [1, 0, 0]

# ❌ INCORRECT
miller_indices = [(1, 0, 0), (1, 1, 0)]
```

### 3. Slab Generation Options

The Ag2O builder includes additional slab generation options not in Ag3PO4:

- `lll_reduce`: Apply LLL reduction to find shorter lattice vectors
- `center_slab`: Center slab in the cell
- `symmetrize`: Symmetrize the slab
- `primitive`: Use primitive cell
- `in_unit_planes`: Express thickness in unit planes instead of Angstroms
- `max_normal_search`: Maximum search for surface normal

These options provide more control over slab generation.

### 4. Surface Energy Output

For binary oxides, surface energies are calculated as:

**γ(Δμ_O)** = Surface energy as a function of oxygen chemical potential

This is different from ternary oxides which give **γ(Δμ_P, Δμ_O)**.

## Complete Example

```python
#!/usr/bin/env python
"""Generate and relax (100) surface of Ag2O (binary oxide)."""

from aiida import load_profile
from teros.core.builders.default_ag2o_builders import get_ag2o_defaults
from teros.core.workgraph import build_core_workgraph_with_map

# Load AiiDA
load_profile(profile='psteros')

# Get defaults for Ag2O (binary oxide)
defaults = get_ag2o_defaults(
    structures_dir="/home/thiagotd/git/PS-TEROS/examples/structures",
    code_label="VASP-VTST-6.4.3@bohr",
    potential_family="PBE",
    miller_indices=[1, 0, 0],  # (100) surface
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    symmetrize=True
)

# Optional: Override specific parameters
defaults['bulk_parameters']['ENCUT'] = 600
defaults['thermodynamics_sampling'] = 200

# Create and submit workgraph
wg = build_core_workgraph_with_map(
    **defaults,
    name="Ag2O_100_binary_oxide"
)

wg.submit(wait=False)

print(f"Submitted Ag2O (binary oxide) calculation")
print(f"WorkGraph PK: {wg.pk}")
print(f"Surface: (100)")
print(f"System type: Binary oxide (Ag2O)")
print(f"Surface energy: γ(Δμ_O)")
```

## Example File

A complete working example is available:

```
examples/default_builders/slabs_autogen_ag2o.py
```

Run it with:

```bash
cd /home/thiagotd/git/worktree/PS-TEROS/default-builders
source ~/envs/psteros/bin/activate
python examples/default_builders/slabs_autogen_ag2o.py
```

## Accessing Results

After the calculation completes:

```python
from aiida import load_node

# Load the workgraph
wg = load_node(PK)  # Replace PK with your workgraph PK

# Get formation enthalpy (for Ag2O)
formation_enthalpy = wg.outputs.formation_enthalpy.get_dict()
print(f"ΔH_f(Ag2O) = {formation_enthalpy['formation_enthalpy']} eV/formula")

# Get generated slabs
slabs = wg.outputs.slab_structures
print(f"Generated terminations: {list(slabs.keys())}")

# Get relaxed slabs
relaxed = wg.outputs.relaxed_slabs
for term_id in relaxed.keys():
    slab = relaxed[term_id]
    print(f"Relaxed {term_id}: {len(slab.sites)} atoms")

# Get surface energies (binary oxide: γ(Δμ_O))
surface_energies = wg.outputs.surface_energies
for term_id in surface_energies.keys():
    data = surface_energies[term_id].get_dict()
    print(f"{term_id}:")
    print(f"  γ(Δμ_O) range: {min(data['gamma_array'])} to {max(data['gamma_array'])} J/m²")
    print(f"  Δμ_O range: {min(data['delta_mu_O'])} to {max(data['delta_mu_O'])} eV")

# Export relaxed slab
atoms = relaxed['term_0'].get_ase()
atoms.write('ag2o_relaxed_100.cif')
```

## Comparison: Binary vs Ternary

| Feature | Ag₂O (Binary) | Ag₃PO₄ (Ternary) |
|---------|---------------|------------------|
| **Elements** | Ag + O | Ag + P + O |
| **Nonmetal** | Dummy (uses Ag) | Physical (P) |
| **Surface Energy** | γ(Δμ_O) | γ(Δμ_P, Δμ_O) |
| **Parameters** | 37 total | 31 total |
| **Slab Options** | 9 options | 3 options |
| **Module** | `default_ag2o_builders` | `default_ag3po4_builders` |

## Validation

Test the builder:

```bash
cd /home/thiagotd/git/worktree/PS-TEROS/default-builders
source ~/envs/psteros/bin/activate

python -c "
from teros.core.builders.default_ag2o_builders import get_ag2o_defaults
d = get_ag2o_defaults(structures_dir='/tmp', code_label='test', potential_family='PBE')
print(f'Parameters: {len(d)}')
print(f'Miller index: {d[\"miller_indices\"]}')
print(f'Bulk: {d[\"bulk_name\"]}')
print(f'Nonmetal (dummy): {d[\"nonmetal_name\"]}')
print(f'Binary oxide: Ready!')
"
```

## Summary

✅ **Module**: `teros.core.builders.default_ag2o_builders`
✅ **Parameters**: 37 (complete set for binary oxides)
✅ **System Type**: Binary oxide (Ag + O only)
✅ **Nonmetal**: Dummy reference (uses Ag parameters)
✅ **Surface Energy**: γ(Δμ_O)
✅ **Slab Generation**: Full control with 9 options
✅ **Example**: `examples/default_builders/slabs_autogen_ag2o.py`

---

**Created**: 2025-01-09
**Module**: `teros.core.builders.default_ag2o_builders`
**System**: Binary Oxide (Ag₂O)
**Status**: ✅ Ready for use
