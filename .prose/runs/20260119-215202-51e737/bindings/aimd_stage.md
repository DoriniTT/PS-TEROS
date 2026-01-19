# AIMD Stage Extraction Complete

## Summary

Successfully extracted AIMD calculation logic from `build_core_workgraph()` (lines 1698-1862) into a new stage module.

## Created File

**File:** `/home/thiagotd/git/PS-TEROS/teros/core/stages/aimd_stage.py`

## Functions Implemented

### 1. `resolve_aimd_parameters()`
Resolves AIMD parameters with fallback chain: aimd -> slab -> bulk.

**Parameters:**
- `aimd_parameters`, `aimd_options`, `aimd_potential_mapping`, `aimd_kpoints_spacing`
- `slab_parameters`, `slab_options`, `slab_potential_mapping`, `slab_kpoints_spacing`
- `bulk_parameters`, `bulk_options`, `bulk_potential_mapping`, `kpoints_spacing`

**Returns:** Tuple of (params, opts, pot_map, kpts)

### 2. `prepare_fixed_atoms()`
Prepares fixed atoms lists for each slab structure.

**Parameters:**
- `input_slabs`: Dict of input slab structures
- `fix_atoms`, `fix_type`, `fix_thickness`, `fix_elements`
- `log`: Logger instance

**Returns:** Dict mapping slab labels to lists of 1-based atom indices

### 3. `determine_initial_slabs_source()`
Determines the source of initial structures for AIMD.

**Priority order:**
1. Relaxed slabs from `relax_slabs_scatter` task
2. Relaxed slabs from `collect_slab_outputs_restart` task
3. Unrelaxed input slabs
4. Generated slabs from `generate_slab_structures` task

**Returns:** Either a dict or task output socket

### 4. `add_aimd_stage()` (Main Function)
Adds the complete AIMD stage to the workgraph.

**Handles:**
1. Determining calculator type (VASP/CP2K)
2. Loading appropriate code
3. Handling fixed atoms if requested
4. Determining initial slab source
5. Creating supercells if requested
6. Creating sequential AIMD stage tasks
7. Wiring outputs between stages

**Parameters:**
- `wg`: WorkGraph instance
- `calculator`: 'vasp' or 'cp2k'
- `aimd_sequence`: List of stage configurations
- `code_label`, `aimd_code_label`
- `aimd_parameters`, `aimd_options`, `aimd_potential_mapping`, `aimd_kpoints_spacing`
- `potential_family`, `clean_workdir`, `max_concurrent_jobs`
- `fix_atoms`, `fix_type`, `fix_thickness`, `fix_elements`, `fix_components`
- `aimd_supercell`, `relax_slabs`, `input_slabs`
- `basis_file`, `pseudo_file` (for CP2K)
- `log`: Logger instance

**Returns:** List of stage task objects

## Logic Preserved

All existing logic from the original implementation was preserved:

1. **Calculator dispatch:** Imports appropriate scatter function based on calculator type
2. **Code loading:** Uses `aimd_code_label` if provided, falls back to `code_label`
3. **Fixed atoms handling:** Pre-computes for input_slabs, supports dynamic calculation for generated slabs
4. **Initial slabs source:** Checks multiple task sources in priority order
5. **Supercell creation:** Adds `create_aimd_supercells` task if `aimd_supercell` is provided
6. **Sequential stages:** Creates AIMD stage tasks with proper output wiring
7. **VASP vs CP2K inputs:** Different input parameters based on calculator type

## Dependencies

```python
from teros.core.aimd_functions import aimd_single_stage_scatter  # VASP
from teros.core.aimd_cp2k import aimd_single_stage_scatter_cp2k  # CP2K
from teros.core.aimd.tasks import create_supercells_scatter
from teros.core.fixed_atoms import get_fixed_atoms_list
from aiida.orm import load_code
```

## Exports Added to `__all__`

```python
__all__ = [
    "resolve_aimd_parameters",
    "prepare_fixed_atoms",
    "determine_initial_slabs_source",
    "add_aimd_stage",
]
```

## Usage Example

```python
from teros.core.stages.aimd_stage import add_aimd_stage, resolve_aimd_parameters

# Resolve parameters
params, opts, pot_map, kpts = resolve_aimd_parameters(...)

# Add AIMD stage to workgraph
stage_tasks = add_aimd_stage(
    wg=wg,
    calculator='vasp',
    aimd_sequence=[
        {'TEBEG': 300, 'NSW': 1000},
        {'TEBEG': 500, 'NSW': 2000},
    ],
    code_label='VASP-6.5.1@cluster',
    aimd_parameters=params,
    aimd_options=opts,
    aimd_potential_mapping=pot_map,
    aimd_kpoints_spacing=kpts,
    ...
)
```

## Notes

- The module follows the existing stage pattern from `restart_handling.py`
- Google-style docstrings with type hints are used throughout
- All imports are conditional/lazy to avoid circular dependencies
- The `add_aimd_stage()` function is the most complex stage function, handling both VASP and CP2K calculators with different parameter requirements
