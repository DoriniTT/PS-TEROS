# Cleavage Energy Calculations

## Overview

The cleavage energy module (`teros.core.cleavage`) calculates the energy required to split a crystal into two complementary surfaces. This is an important property for understanding crystal stability and surface reactivity.

## Theory

The cleavage energy for a pair of complementary slabs is calculated using:

```
Ec(i,j) = 1/(2A) * (E_i^slab + E_j^slab - n*E_bulk)
```

Where:
- `E_i^slab`, `E_j^slab`: Total energies of complementary slabs (eV)
- `A`: Surface area of the slab (Ų)
- `n`: Number of formula units in the combined slabs
- Factor `1/(2A)`: Accounts for two surfaces created by cleavage

The result is given in both eV/Ų and J/m² (1 eV/Ų = 16.0217663 J/m²).

## Complementary Slab Pairs

Pymatgen's SlabGenerator produces slabs in a consistent complementary ordering:

### Even number of terminations (e.g., 4)
- Termination 0 is complementary to termination 3
- Termination 1 is complementary to termination 2
- Pairs: (0,3), (1,2)

### Odd number of terminations (e.g., 5)
- Termination 0 is complementary to termination 4
- Termination 1 is complementary to termination 3
- Termination 2 is complementary to itself
- Pairs: (0,4), (1,3), (2,2)

## Usage in WorkGraph

### Enable cleavage energy calculations

Add the following parameter when building your WorkGraph:

```python
wg = build_core_workgraph_with_map(
    # ... other parameters ...
    relax_slabs=True,  # Required: must relax slabs first
    compute_cleavage=True,  # Enable cleavage energy calculations
)
```

**Note**: Cleavage energy calculations require `relax_slabs=True` because they need the relaxed slab structures and their energies.

### Access results

After the WorkGraph completes:

```python
from aiida import load_node

wg = load_node(PK)

# Access cleavage energies
cleavages = wg.outputs.cleavage_energies

# List all pairs
print(list(cleavages.keys()))  # e.g., ['pair_0_3', 'pair_1_2']

# Get data for a specific pair
pair_data = cleavages['pair_0_3'].get_dict()

# Cleavage energy values
ec_ev = pair_data['cleavage_energy_eV_A2']  # eV/Ų
ec_jm2 = pair_data['cleavage_energy_J_m2']  # J/m²

# Additional information
area = pair_data['area_A2']  # Surface area
n_fu = pair_data['n_formula_units']  # Formula units
term_i = pair_data['term_i']  # First termination index
term_j = pair_data['term_j']  # Second termination index
```

## Output Data Structure

Each cleavage energy result (per pair) contains:

- `cleavage_energy_eV_A2`: Cleavage energy in eV/Ų
- `cleavage_energy_J_m2`: Cleavage energy in J/m²
- `term_i`, `term_j`: Indices of the complementary terminations
- `pair_label`: String label for the pair (e.g., "term_0_term_3")
- `n_formula_units`: Number of formula units in combined slabs
- `area_A2`: Surface area in ų
- `E_slab_i_eV`, `E_slab_j_eV`: Total energies of individual slabs
- `E_bulk_eV`: Total energy of bulk structure
- `bulk_stoichiometry`: Reduced stoichiometry per formula unit
- `slab_i_composition`, `slab_j_composition`: Atomic composition of each slab
- `combined_composition`: Combined atomic composition

## Works for Both Binary and Ternary Oxides

The module automatically handles:

### Binary oxides (e.g., Ag₂O)
- Bulk formula: Ag₂O
- Reduced stoichiometry: 2 Ag + 1 O per formula unit
- Calculates n from combined slab composition

### Ternary oxides (e.g., Ag₃PO₄)
- Bulk formula: Ag₃PO₄
- Reduced stoichiometry: 3 Ag + 1 P + 4 O per formula unit
- Calculates n from combined slab composition

## Example

See `examples/cleavage/slabs_relax_ag2o_cleavage.py` for a complete working example with Ag₂O.

## Reference

Yang et al., "Cleavage energy for symmetric surfaces" formula:
```
Ec(i,j) = 1/(2A) * (E_i^slab + E_j^slab - n*E_bulk)
```

This approach is particularly useful for understanding surface reconstruction energies and comparing stability of different surface orientations.
