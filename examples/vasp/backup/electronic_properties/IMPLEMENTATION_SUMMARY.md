# Implementation Summary: B-based Formulation for Ternary Oxides

## Overview

The `calculate_surface_energy_ternary` function in `teros/core/thermodynamics.py` has been extended to compute surface free energy with respect to **both** metal elements in a ternary oxide system.

## What Was Added

### Dual Formulation Support

For a ternary oxide **A_x B_y O_z** (e.g., Ag₃PO₄), the function now returns:

1. **A-based formulation**: γ(Δμ_A, Δμ_O)
   - Element A as independent variable (e.g., Ag)
   - Element B as reference (eliminated via bulk equilibrium) (e.g., P)
   
2. **B-based formulation**: γ(Δμ_B, Δμ_O)  [NEW]
   - Element B as independent variable (e.g., P)
   - Element A as reference (eliminated via bulk equilibrium) (e.g., Ag)

### Output Structure

The function returns a dictionary with three main sections:

```python
{
    'A_based': {
        'phi': float,                    # Reference surface energy
        'Gamma_A': float,                # Surface excess of A
        'Gamma_O': float,                # Surface excess of O
        'delta_mu_A_range': list,        # Chemical potential range for A
        'delta_mu_O_range': list,        # Chemical potential range for O
        'gamma_grid': 2D list,           # Surface energy γ(Δμ_A, Δμ_O)
        'gamma_at_muA_zero': list,       # 1D slice at Δμ_A=0
        'gamma_at_muO_zero': list,       # 1D slice at Δμ_O=0
        'element_A_independent': str,    # Name of element A
        'element_B_reference': str,      # Name of element B
    },
    
    'B_based': {
        'phi': float,                    # Reference surface energy (B-based)
        'Gamma_B': float,                # Surface excess of B
        'Gamma_O': float,                # Surface excess of O (with A as ref)
        'delta_mu_B_range': list,        # Chemical potential range for B
        'delta_mu_O_range': list,        # Chemical potential range for O
        'gamma_grid': 2D list,           # Surface energy γ(Δμ_B, Δμ_O)
        'gamma_at_muB_zero': list,       # 1D slice at Δμ_B=0
        'gamma_at_muO_zero': list,       # 1D slice at Δμ_O=0
        'element_B_independent': str,    # Name of element B
        'element_A_reference': str,      # Name of element A
    },
    
    # Common information
    'oxide_type': 'ternary',
    'area_A2': float,
    'bulk_stoichiometry_AxByOz': dict,
    'slab_atom_counts': dict,
    'reference_energies_per_atom': dict,
    'E_slab_eV': float,
    'E_bulk_fu_eV': float,
    'formation_enthalpy_eV': float,
    
    # Legacy keys (for backward compatibility)
    'phi': float,
    'Gamma_M_vs_Nref': float,
    'Gamma_O_vs_Nref': float,
    'delta_mu_M_range': list,
    'delta_mu_O_range': list,
    'gamma_grid': 2D list,
    'gamma_at_muM_zero': list,
    'gamma_at_muO_zero': list,
    'element_M_independent': str,
    'element_N_reference': str,
}
```

## Mathematical Framework

### A-based Formulation (Original)

For ternary oxide A_x B_y O_z, eliminate μ_B using bulk equilibrium:

```
μ_B = (E_bulk - x·μ_A - z·μ_O) / y
```

Surface energy becomes:
```
γ(Δμ_A, Δμ_O) = φ_A - Γ_A·Δμ_A - Γ_O·Δμ_O
```

where:
- Γ_A = (N_A - x·N_B/y) / (2A)
- Γ_O = (N_O - z·N_B/y) / (2A)

### B-based Formulation (New)

Eliminate μ_A instead of μ_B:

```
μ_A = (E_bulk - y·μ_B - z·μ_O) / x
```

Surface energy becomes:
```
γ(Δμ_B, Δμ_O) = φ_B - Γ_B·Δμ_B - Γ_O·Δμ_O
```

where:
- Γ_B = (y·N_A/x - N_B) / (2A)
- Γ_O = (z·N_A/x - N_O) / (2A)

## Use Cases

### For Ag₃PO₄ System

**A-based (original):** Analyze surface stability as a function of silver chemical potential
- Useful for understanding Ag-rich vs Ag-poor environments
- Γ_Ag indicates silver excess/deficiency at surface
- Use when interested in silver deposition/dissolution

**B-based (new):** Analyze surface stability as a function of phosphorus chemical potential
- Useful for understanding P-rich vs P-poor environments
- Γ_P indicates phosphorus excess/deficiency at surface
- Use when interested in phosphate chemistry effects

## Backward Compatibility

The function maintains full backward compatibility by including legacy keys at the top level of the output dictionary. Existing code that accesses `result['gamma_grid']` or `result['phi']` will continue to work unchanged, receiving the A-based formulation results.

## Theoretical Documentation

See `derivation_B_element.tex` for the complete mathematical derivation of the B-based formulation, including:
- Detailed step-by-step derivation
- Chemical potential constraints
- Comparison with A-based formulation
- Application to Ag₃PO₄

## Files Modified

1. **teros/core/thermodynamics.py**
   - Modified `calculate_surface_energy_ternary()` function
   - Added ~80 lines of code for B-based computation
   - Maintained backward compatibility

2. **derivation_B_element.tex** (new)
   - Complete LaTeX documentation of the derivation
   - Can be compiled to PDF with: `pdflatex derivation_B_element.tex`

## Testing

The implementation has been tested with mock Ag₃PO₄ structures and successfully produces:
- Both A-based and B-based results
- Correct surface excesses (Γ_A, Γ_B, Γ_O)
- Valid chemical potential ranges
- 2D gamma grids for plotting phase diagrams

## Next Steps

To use this in your workflows:

1. **Access A-based results** (Ag as independent):
   ```python
   gamma_Ag_O = result['A_based']['gamma_grid']
   delta_mu_Ag = result['A_based']['delta_mu_A_range']
   ```

2. **Access B-based results** (P as independent):
   ```python
   gamma_P_O = result['B_based']['gamma_grid']
   delta_mu_P = result['B_based']['delta_mu_B_range']
   ```

3. **Plot both phase diagrams** to understand surface stability from different perspectives

4. **Compare surface excesses** to understand which element is more important for surface structure
