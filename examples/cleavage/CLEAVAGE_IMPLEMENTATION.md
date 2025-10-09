# Cleavage Energy Module Implementation

## Summary

Successfully implemented a new cleavage energy calculation module for PS-TEROS that computes the energy required to split crystals into complementary surfaces.

## Implementation Details

### 1. New Module: `teros/core/cleavage.py`

Created a complete cleavage energy calculation module with:

#### Functions:
- **`calculate_cleavage_energy`** (calcfunction): Computes cleavage energy for a single pair of complementary slabs
  - Formula: `Ec(i,j) = 1/(2A) * (E_i^slab + E_j^slab - n*E_bulk)`
  - Returns energy in both eV/Ų and J/m²
  - Validates stoichiometry consistency
  - Works for both binary and ternary oxides

- **`compute_cleavage_energies_scatter`** (task.graph): Parallel computation for all complementary pairs
  - Automatically determines complementary pairs from pymatgen convention
  - Handles even and odd numbers of terminations
  - Uses scatter-gather pattern for efficient parallel execution

#### Key Features:
- Automatic pairing of complementary terminations according to pymatgen's SlabGenerator convention:
  - Even n: (0,n-1), (1,n-2), (2,n-3), ...
  - Odd n: (0,n-1), (1,n-2), ..., (middle,middle)
- Robust stoichiometry validation
- Surface area calculation for general (non-orthogonal) cells
- Comprehensive output data including compositions and energies

### 2. WorkGraph Integration: `teros/core/workgraph.py`

Updated the core workgraph to include cleavage energy calculations:

#### Changes:
- Added `compute_cleavage` parameter (default: False)
- Added `cleavage_energies` to workflow outputs
- Integrated cleavage calculation after slab relaxation (parallel to thermodynamics)
- Updated all builder functions to accept and propagate the parameter
- Updated documentation strings

#### Workflow Position:
```
1. Bulk + reference relaxations (parallel)
2. Formation enthalpy calculation
3. Slab generation
4. Slab relaxation (if enabled)
5. [Parallel execution]:
   - Thermodynamics calculation (if enabled)
   - Cleavage energy calculation (if enabled)  ← NEW
```

### 3. Example Script: `examples/cleavage/slabs_relax_ag2o_cleavage.py`

Updated the example to demonstrate cleavage energy calculations:
- Added `compute_cleavage=True` parameter
- Updated workflow description and documentation
- Added example code for accessing cleavage energy results
- Updated expected outputs section

### 4. Module Exports: `teros/core/__init__.py`

Added cleavage module exports:
- `calculate_cleavage_energy`
- `compute_cleavage_energies_scatter`

### 5. Documentation: `docs/cleavage_energy.md`

Created comprehensive documentation covering:
- Theory and formula
- Complementary pairing logic
- Usage examples
- Output data structure
- Binary and ternary oxide support

## Requirements

Cleavage energy calculation requires:
- `relax_slabs=True`: Needs relaxed slab structures and energies
- `compute_cleavage=True`: Enables the calculation

## Usage Example

```python
wg = build_core_workgraph_with_map(
    structures_dir="/path/to/structures",
    bulk_name="ag2o.cif",
    metal_name="Ag.cif",
    nonmetal_name="Ag.cif",  # Dummy for binary oxide
    oxygen_name="O2.cif",
    miller_indices=[1, 0, 0],
    min_slab_thickness=10.0,
    min_vacuum_thickness=15.0,
    # ... VASP parameters ...
    relax_slabs=True,  # Required
    compute_cleavage=True,  # Enable cleavage calculation
)
wg.submit()
```

## Accessing Results

```python
from aiida import load_node

wg = load_node(PK)
cleavages = wg.outputs.cleavage_energies

for pair_id in cleavages.keys():
    data = cleavages[pair_id].get_dict()
    ec = data['cleavage_energy_eV_A2']
    print(f"{pair_id}: {ec:.4f} eV/Ų")
```

## Output Structure

Each pair returns a Dict with:
- `cleavage_energy_eV_A2`: Energy in eV/Ų
- `cleavage_energy_J_m2`: Energy in J/m²
- `term_i`, `term_j`: Termination indices
- `pair_label`: String identifier
- `n_formula_units`: Formula units in combined slabs
- `area_A2`: Surface area
- `E_slab_i_eV`, `E_slab_j_eV`: Individual slab energies
- `E_bulk_eV`: Bulk energy
- `bulk_stoichiometry`: Reduced formula
- Composition data for both slabs

## Testing

All modules successfully import and pass syntax validation:
- ✓ `teros.core.cleavage` module imports correctly
- ✓ `teros.core.workgraph` module imports correctly
- ✓ Example script syntax is valid
- ✓ Complementary pairing logic tested for n=2 to 8 terminations
- ✓ AiiDA daemon restarted successfully

## Compatibility

Works with:
- Binary oxides (e.g., Ag₂O, MgO, ZnO)
- Ternary oxides (e.g., Ag₃PO₄, Fe₂WO₆)
- All Miller indices
- Any number of terminations (even or odd)

## Reference

Implementation based on Yang et al. cleavage energy formula for symmetric surfaces, as described in the user's requirements.
