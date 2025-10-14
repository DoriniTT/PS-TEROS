# Using Pre-Generated Slab Structures in PS-TEROS

This guide explains how to use pre-generated slab structures as input to PS-TEROS, giving you full flexibility to provide your own slabs instead of relying on automatic slab generation.

## Overview

PS-TEROS now supports two modes for working with slab structures:

1. **Automatic Generation Mode** (default): PS-TEROS generates slabs from the relaxed bulk structure using Pymatgen's SlabGenerator
2. **User-Provided Mode** (new): Users provide their own pre-generated slab structures, and PS-TEROS only performs relaxation and energy calculations

## Benefits of User-Provided Slabs

- **Full control**: Use your own slab generation tools or manual construction
- **Flexibility**: Provide slabs with specific surface reconstructions, adsorbates, or defects
- **Reproducibility**: Use exact structures from literature or previous calculations
- **Efficiency**: Skip slab generation if you already have the structures you need
- **Custom terminations**: Work with non-standard or manually edited surface terminations

## How to Use User-Provided Slabs

### 1. Prepare Your Slab Structures

Create your slab structure files in any format supported by ASE (CIF, POSCAR, xyz, etc.). Place them in a directory, for example:

```
examples/slabs/input_structures/
├── slab_term_0.cif
├── slab_term_1.cif
└── slab_term_2.cif
```

**Important considerations for slab structures:**
- Slabs should have sufficient vacuum spacing (typically 10-15 Å)
- The c-axis should be perpendicular to the surface
- Atoms should be centered in the c direction if needed
- Cell parameters should be appropriate for your system

### 2. Load Structures in Your Script

```python
from aiida import orm
from ase.io import read

# Create a dictionary of slab structures
input_slabs = {}

slab_files = ["slab_term_0.cif", "slab_term_1.cif", "slab_term_2.cif"]

for idx, slab_file in enumerate(slab_files):
    slab_path = f"path/to/slabs/{slab_file}"
    atoms = read(slab_path)
    structure = orm.StructureData(ase=atoms)
    structure.store()  # IMPORTANT: Must store before passing!
    input_slabs[f"term_{idx}"] = structure
```

**Key points:**
- Dictionary keys should follow the pattern `"term_0"`, `"term_1"`, etc.
- Each value should be an AiiDA `StructureData` node
- **CRITICAL**: Structures MUST be stored (`.store()`) before adding to dict
- You can load any number of slabs

### 3. Call the WorkGraph Builder

```python
from teros.core.workgraph import build_core_workgraph_with_map

wg = build_core_workgraph_with_map(
    structures_dir=structures_dir,
    bulk_name="ag3po4.cif",
    metal_name="Ag.cif",
    nonmetal_name="P.cif",
    oxygen_name="O2.cif",
    code_label=code_label,
    potential_family=potential_family,
    # ... other parameters ...
    
    # Provide your pre-generated slabs
    input_slabs=input_slabs,
    
    # Enable slab relaxation
    relax_slabs=True,
    slab_parameters=slab_parameters,
    slab_options=slab_options,
    slab_potential_mapping=slab_potential_mapping,
    slab_kpoints_spacing=slab_kpoints_spacing,
    
    name="MyProject_InputSlabs",
)
```

**When using `input_slabs`:**
- `miller_indices`, `min_slab_thickness`, and `min_vacuum_thickness` are **not required**
- Slab generation is automatically skipped
- Your provided structures are used directly

### 4. Submit and Monitor

```python
wg.submit(wait=False)
print(f"WorkGraph PK: {wg.pk}")

# Monitor progress
# verdi process show <PK>
```

## Complete Example

See `examples/slabs/slabs_input_relax.py` for a complete working example that demonstrates:
- Loading pre-generated slab structures from files
- Setting up relaxation parameters
- Creating and submitting the workflow
- Accessing results

## Parameters Comparison

### Automatic Generation Mode

```python
wg = build_core_workgraph_with_map(
    # ... bulk and reference parameters ...
    
    # Required for automatic generation:
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    
    # Optional generation parameters:
    lll_reduce=True,
    center_slab=True,
    symmetrize=True,
    primitive=True,
    
    relax_slabs=True,
    # ... relaxation parameters ...
)
```

### User-Provided Mode

```python
wg = build_core_workgraph_with_map(
    # ... bulk and reference parameters ...
    
    # Provide pre-generated slabs:
    input_slabs=my_slab_dict,
    
    # Generation parameters not needed (will be ignored)
    # miller_indices, min_slab_thickness, etc. not required
    
    relax_slabs=True,
    # ... relaxation parameters ...
)
```

## Outputs

Both modes produce the same output structure:

```python
from aiida import load_node

wg = load_node(PK)

# Input/generated slabs
slabs = wg.outputs.slab_structures
# Keys: 'term_0', 'term_1', ...

# Relaxed slabs (if relax_slabs=True)
relaxed = wg.outputs.relaxed_slabs
# Keys: 'term_0', 'term_1', ...

# Energies (if relax_slabs=True)
energies = wg.outputs.slab_energies
# Keys: 'term_0', 'term_1', ...

# Access specific slab
term_0_energy = energies['term_0'].value
term_0_structure = relaxed['term_0']
```

## Tips and Best Practices

1. **Consistent naming**: Use sequential naming (`term_0`, `term_1`, ...) for better organization
2. **Validate structures**: Check your slab structures before running calculations
3. **Cell parameters**: Ensure cell parameters are appropriate for periodic DFT calculations
4. **Vacuum spacing**: Verify adequate vacuum to avoid slab-slab interactions
5. **Initial geometry**: Starting from good initial structures can improve convergence
6. **Mixed approach**: You can generate some slabs automatically and others manually

## Troubleshooting

**Q: Can I mix automatically-generated and user-provided slabs?**
A: No, you must choose one mode or the other per workflow run. However, you can run separate workflows.

**Q: What file formats are supported?**
A: Any format that ASE's `read()` function supports: CIF, POSCAR, CONTCAR, xyz, etc.

**Q: Do I need to relax my input slabs?**
A: No, you can set `relax_slabs=False` to skip relaxation and use your structures as-is.

**Q: Can I still calculate surface energies with user-provided slabs?**
A: Yes! Set `compute_thermodynamics=True` and `relax_slabs=True` as usual.

## Related Examples

- `examples/backup/slabs/slabs_relax.py` - Original automatic generation example
- `examples/slabs/slabs_input_relax.py` - User-provided slabs example

## Code Changes

The feature was implemented by:
1. Adding `input_slabs` parameter to `core_workgraph()` in `teros/core/workgraph.py`
2. Making slab generation parameters optional when `input_slabs` is provided
3. Conditional logic to use provided slabs or generate them automatically
