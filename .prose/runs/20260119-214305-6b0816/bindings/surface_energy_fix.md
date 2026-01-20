# surface_energy_fix

kind: let

source:
```prose
surface_energy_fix = session: fixer
  prompt: """Update teros/core/surface_hydroxylation/surface_energy.py to use modern task decorators."""
```

---

## Summary of Changes

Updated `teros/core/surface_hydroxylation/surface_energy.py` to use modern AiiDA-WorkGraph task decorators.

### Changes Made

1. **Import statement** (line 9):
   - OLD: `from aiida.engine import calcfunction`
   - NEW: `from aiida_workgraph import task`

2. **Decorator updates** - Changed 8 functions from `@calcfunction` to `@task.calcfunction`:
   - `calc_delta_g_general_reaction1` (line 243)
   - `calc_delta_g_general_reaction2` (line 275)
   - `calc_delta_g_general_reaction3` (line 311)
   - `calc_delta_g_reaction1` (line 346)
   - `calc_delta_g_reaction2` (line 396)
   - `calc_delta_g_reaction3` (line 436)
   - `calc_gamma_s` (line 541)
   - `calc_gamma` (line 577)

### Preserved

- All docstrings remain unchanged
- All function logic remains unchanged
- All helper functions (non-decorated) remain unchanged:
  - `_get_formula_dict`
  - `analyze_composition`
  - `analyze_composition_general`
  - `_calc_delta_g_general_reaction1_internal`
  - `_calc_delta_g_general_reaction2_internal`
  - `_calc_delta_g_general_reaction3_internal`
  - `select_reaction_function`
  - `get_surface_area`
  - `calculate_surface_energies_general`
  - `calculate_surface_energies`
