# Adsorption Stage Extraction

## Status: COMPLETED

## File Created
- `/home/thiagotd/git/PS-TEROS/teros/core/stages/adsorption_stage.py`

## Functions Implemented

### 1. `resolve_adsorption_parameters()`
Resolves adsorption parameters with fallback chain (adsorption -> slab -> bulk).

**Parameters:**
- `adsorption_parameters`, `adsorption_options`, `adsorption_potential_mapping`, `adsorption_kpoints_spacing`
- `slab_parameters`, `slab_options`, `slab_potential_mapping`, `slab_kpoints_spacing`
- `bulk_parameters`, `bulk_options`, `bulk_potential_mapping`, `kpoints_spacing`

**Returns:** Tuple of (resolved_params, resolved_opts, resolved_pot_map, resolved_kpts)

### 2. `extract_builder_params()`
Extracts INCAR parameters from builder_inputs for backward compatibility.

**Parameters:**
- `relax_builder_inputs`: New-style builder inputs for relaxation
- `scf_builder_inputs`: New-style builder inputs for SCF
- `fallback_params`: Fallback parameters (old-style)

**Returns:** Tuple of (relax_params, scf_params)

### 3. `add_adsorption_energy_task()`
Adds the compute_adsorption_energies_scatter task and connects outputs.

**Parameters:**
- `wg`: WorkGraph instance
- `code`: Loaded AiiDA Code
- All adsorption-related parameters

**Connected outputs:**
- `relaxed_complete_structures`
- `separated_structures`
- `substrate_energies`
- `molecule_energies`
- `complete_energies`
- `adsorption_energies`

### 4. `add_adsorption_energy_stage()`
Main orchestrator for Stage 10. Performs:
1. Load VASP code
2. Resolve parameters with fallback chain
3. Extract relax/SCF parameters from builder_inputs
4. Add compute_adsorption_energies_scatter task
5. Connect all outputs to wg.outputs

**Parameters:**
- `wg`: WorkGraph instance
- `code_label`: VASP code label
- `adsorption_structures`: Dict of structures
- `adsorption_formulas`: Dict of adsorbate formulas
- `potential_family`: POTCAR family
- All adsorption-specific parameters (optional, with fallbacks)
- Atom fixing parameters
- Fallback parameters (slab/bulk)
- `logger`: Optional logger

## Source Lines Extracted
Lines 1864-1965 from `build_core_workgraph()` in `workgraph.py`

## Patterns Followed
- Google-style docstrings with Args, Returns, Example sections
- Type hints using `Optional`, `Dict`, `Any`, `List` from typing
- TYPE_CHECKING guard for AiiDA imports
- Module-level logger fallback
- Consistent with existing stage modules (restart_handling.py, preset_resolution.py)
