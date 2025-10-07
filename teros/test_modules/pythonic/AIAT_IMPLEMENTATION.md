# Ab Initio Atomistic Thermodynamics Implementation

This document describes the integration of ab initio atomistic thermodynamics (AIAT) calculations into the pythonic scatter-gather workflow for ternary oxide surfaces.

## Implementation Summary

✅ **Created**: `aiat_ternary.py` - Complete AIAT module
✅ **Updated**: `workgraph.py` - Added thermodynamics step
✅ **Updated**: `slabs_relax.py` - Added `--with-thermodynamics` flag

## Files

### aiat_ternary.py

New module containing:

1. **`calculate_surface_energy_ternary`** (`@task.calcfunction`)
   - Computes γ(Δμ_M, Δμ_O) for a single slab
   - Returns surface energy data including φ, Γ_M, Γ_O
   - Calculates energy grid over chemical potential range

2. **`compute_surface_energies_scatter`** (`@task.graph`)
   - Scatter-gather pattern for computing all slab surface energies in parallel
   - Iterates over relaxed slabs and their energies
   - Returns dynamic namespace of surface energy data

3. **Mock data generators** (for testing):
   - `create_mock_reference_energies()` 
   - `create_mock_bulk_energy()`
   - `create_mock_formation_enthalpy()`

### workgraph.py Updates

**`slab_relaxation_scatter_gather`** now accepts:
- `compute_thermodynamics`: bool flag
- `bulk_energy`: orm.Float (bulk total energy)
- `reference_energies`: orm.Dict (element reference energies)
- `formation_enthalpy`: orm.Float (formation enthalpy)
- `sampling`: int (grid resolution)

Returns additional output:
- `surface_energies`: Dynamic namespace of AIAT data

### slabs_relax.py Updates

New command-line arguments:
- `--with-thermodynamics`: Enable AIAT calculations
- `--sampling`: Grid resolution (default: 100)

## Usage

### Mock Mode with Thermodynamics

```bash
python slabs_relax.py --mock --with-thermodynamics --mock-count 2
```

Output includes thermodynamics results for each mock value.

### Full VASP Workflow with Thermodynamics

```bash
python slabs_relax.py --with-thermodynamics --sampling 100
```

Requires:
- Bulk structure and energy
- Reference energies for elements
- Formation enthalpy

Currently uses mock data for these inputs (see slabs_relax.py lines ~165-175).

## Thermodynamics Theory

For ternary oxide M-N-O system:

**Surface energy**: γ(Δμ_M, Δμ_O) = φ - Γ_M·Δμ_M - Γ_O·Δμ_O

Where:
- **φ**: Reference surface energy at bulk equilibrium
- **Γ_M, Γ_O**: Surface excess of M and O (relative to reference metal N)
- **Δμ_M, Δμ_O**: Chemical potential deviations from reference state

**Surface excess**:
- Γ_M = (N_M - x_M·N_N/y_N) / (2·A)
- Γ_O = (N_O - z_O·N_N/y_N) / (2·A)

Where M_x N_y O_z is the bulk stoichiometry.

## Integration Points

The thermodynamics step runs AFTER slab relaxations complete:

```
1. Generate slabs
2. Relax slabs in parallel (VASP)
3. Extract energies
4. → Compute surface energies in parallel (AIAT)
5. Return all results
```

Each slab's surface energy calculation is independent → parallel execution.

## Outputs

Each slab gets a comprehensive thermodynamics dictionary:

```python
{
    'phi': float,  # Reference surface energy
    'Gamma_M_vs_Nref': float,  # Surface excess M
    'Gamma_O_vs_Nref': float,  # Surface excess O
    'gamma_values_grid': dict,  # γ(Δμ_M, Δμ_O) grid
    'gamma_values_fixed_muM_zero': dict,  # γ with Δμ_M=0
    'area_A2': float,  # Surface area
    'element_M_independent': str,  # M element
    'element_N_reference': str,  # N element
    'bulk_stoichiometry_MxNyOz': dict,
    'slab_atom_counts': dict,
    'reference_energies_per_atom': dict,
    'E_slab_eV': float,
    'E_bulk_fu_eV': float,
    'formation_enthalpy_eV': float,
    'delta_mu_M_range': list,
    'delta_mu_O_range': list,
}
```

## Example Results

For Ag₃PO₄ (100) surface with 4 terminations:

```
verdi process show <PK>

Outputs:
  surface_energies
    slab_00    <PK>  Dict
    slab_01    <PK>  Dict
    slab_02    <PK>  Dict
    slab_03    <PK>  Dict
```

Each Dict contains the full thermodynamics data above.

## Next Steps

To use with real data:

1. **Bulk calculation**: Run VASP on bulk structure, extract total energy
2. **Reference calculations**: Run VASP on elemental references (Ag metal, P, O₂)
3. **Formation enthalpy**: Calculate from bulk and references
4. **Update slabs_relax.py**: Replace mock data generators with real AiiDA nodes

## Testing

Mock mode verified working:
```bash
cd teros/test_modules/pythonic
python slabs_relax.py --mock --with-thermodynamics
# ✅ Generates 2 mock values
# ✅ Shifts them in parallel
# ✅ Computes mock thermodynamics in parallel
# ✅ Returns all results with provenance
```

