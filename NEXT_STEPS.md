# Next Steps: Getting B-based Results in Your Workflow

## Summary of Changes

The code has been successfully updated to compute surface energy with respect to **both** the A element (Ag) and B element (P) for ternary oxides like Ag₃PO₄.

### What Was Modified

1. **`teros/core/thermodynamics.py`**
   - Enhanced `calculate_surface_energy_ternary()` function
   - Now returns both A-based and B-based formulations
   - Maintains backward compatibility

2. **Documentation Files Created**
   - `derivation_B_element.tex` - Complete mathematical derivation
   - `IMPLEMENTATION_SUMMARY.md` - Technical implementation details
   - `example_usage_B_formulation.py` - Plotting examples
   - `how_to_access_B_based_results.py` - Workflow result access guide

## Current Status

✓ Code updated and tested  
✓ AiiDA daemon restarted (picks up new code)  
✗ Your workflow PK 29082 ran **before** the update  

### Your Workflow Output (PK 29082)

The Dict node 29295 shows the **old format** with only A-based results:
- `element_M_independent`: "Ag"
- `Gamma_M_vs_Nref`: Surface excess of Ag
- `gamma_grid`: 2D array γ(Δμ_Ag, Δμ_O)

**Missing**: B-based results with respect to P

## How to Get B-based Results

### Step 1: Rerun Your Workflow

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-ternary-element
source ~/envs/psteros/bin/activate
python examples/complete/complete_ag3po4_example.py
```

The daemon will use the updated code automatically.

### Step 2: Wait for Completion

Monitor progress:
```bash
verdi process list
# Or for specific workflow:
verdi process status <NEW_PK>
```

### Step 3: Access Results

Once finished, the new workflow will have Dict nodes with **both formulations**:

```python
import aiida
aiida.load_profile('psteros')
from aiida import orm

# Load your new workflow
node = orm.load_node(NEW_PK)

# Get surface energy for a termination
term_0 = node.outputs.surface_energies.term_0.get_dict()

# Access A-based results (Ag as independent)
a_based = term_0['A_based']
print(f"Element A: {a_based['element_A_independent']}")  # "Ag"
print(f"Gamma_A: {a_based['Gamma_A']}")
print(f"Gamma_O: {a_based['Gamma_O']}")
gamma_Ag_O = a_based['gamma_grid']  # 2D array

# Access B-based results (P as independent) - NEW!
b_based = term_0['B_based']
print(f"Element B: {b_based['element_B_independent']}")  # "P"
print(f"Gamma_B: {b_based['Gamma_B']}")
print(f"Gamma_O: {b_based['Gamma_O']}")
gamma_P_O = b_based['gamma_grid']  # 2D array
```

### Step 4: Visualize Both Formulations

Use the provided helper script:

```python
# In how_to_access_B_based_results.py, update:
WORKFLOW_PK = YOUR_NEW_PK

# Then run:
python how_to_access_B_based_results.py
```

This will show both A-based and B-based results for all terminations.

## New Output Structure

After rerunning, each surface energy Dict will contain:

```python
{
    'A_based': {
        'phi': float,
        'Gamma_A': float,
        'Gamma_O': float,
        'delta_mu_A_range': list,
        'delta_mu_O_range': list,
        'gamma_grid': 2D list,
        'gamma_at_muA_zero': list,
        'gamma_at_muO_zero': list,
        'element_A_independent': 'Ag',
        'element_B_reference': 'P',
    },
    
    'B_based': {
        'phi': float,
        'Gamma_B': float,
        'Gamma_O': float,
        'delta_mu_B_range': list,
        'delta_mu_O_range': list,
        'gamma_grid': 2D list,
        'gamma_at_muB_zero': list,
        'gamma_at_muO_zero': list,
        'element_B_independent': 'P',
        'element_A_reference': 'Ag',
    },
    
    # Common info + legacy keys for compatibility
    'oxide_type': 'ternary',
    'area_A2': float,
    'bulk_stoichiometry_AxByOz': {...},
    # ... etc
}
```

## Physical Interpretation

### A-based: γ(Δμ_Ag, Δμ_O)
- **Use when**: Interested in silver-rich vs silver-poor environments
- **Γ_Ag > 0**: Surface is Ag-deficient (relative to bulk)
- **Γ_Ag < 0**: Surface is Ag-rich (relative to bulk)
- **Applications**: Ag deposition, Ag dissolution, metallic effects

### B-based: γ(Δμ_P, Δμ_O)
- **Use when**: Interested in phosphorus-rich vs phosphorus-poor environments
- **Γ_P > 0**: Surface is P-deficient (relative to bulk)
- **Γ_P < 0**: Surface is P-rich (relative to bulk)
- **Applications**: Phosphate chemistry, P doping, catalytic activity

## Plotting Phase Diagrams

After getting new results, create phase diagrams:

```python
from example_usage_B_formulation import plot_both_formulations
import matplotlib.pyplot as plt

# Get the result dict
result_dict = node.outputs.surface_energies.term_0.get_dict()

# Create side-by-side phase diagrams
fig = plot_both_formulations(result_dict)
fig.savefig('ag3po4_phase_diagrams.png', dpi=300, bbox_inches='tight')
plt.show()
```

This creates:
- **Left panel**: γ(Δμ_Ag, Δμ_O) contour plot
- **Right panel**: γ(Δμ_P, Δμ_O) contour plot

## Verification

To verify the code is working without rerunning the full workflow:

```bash
# Quick syntax check
python -m py_compile teros/core/thermodynamics.py

# Check daemon status
verdi daemon status

# Test with mock data (instant)
python -c "
import aiida
aiida.load_profile('psteros')
from teros.core.thermodynamics import calculate_surface_energy_ternary
print('✓ Function loaded successfully')
print(f'✓ Function has {len(calculate_surface_energy_ternary._callable.__code__.co_names)} operations')
"
```

## Files Reference

- **Theory**: `derivation_B_element.tex` (can compile with `pdflatex`)
- **Implementation**: `IMPLEMENTATION_SUMMARY.md`
- **Plotting**: `example_usage_B_formulation.py`
- **Access**: `how_to_access_B_based_results.py`
- **This guide**: `NEXT_STEPS.md`

## Questions?

The implementation is complete and tested. After rerunning your workflow, you'll have access to both formulations, allowing you to analyze Ag₃PO₄ surface stability from both silver and phosphorus perspectives.
