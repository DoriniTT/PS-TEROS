# preset_stage

kind: let

source:
```prose
let preset_stage = session: implementer
  context: analysis, stages_init
  prompt: """Extract the workflow preset resolution logic into a stage function..."""
```

---

## Summary

Created the `teros/core/stages/preset_resolution.py` module implementing Stage 1 of the `build_core_workgraph()` decomposition.

### Function Implemented

**`resolve_workflow_preset()`** - Resolves workflow preset and applies user overrides.

This function performs the complete preset resolution process:

1. **Check old-style API**: Emit deprecation warning if user specifies flags without an explicit `workflow_preset` parameter.

2. **Resolve preset**: Load the preset configuration and apply any user overrides for individual flags.

3. **Validate preset inputs**: Check that all required parameters for the selected preset are provided.

4. **Validate flag dependencies**: Check that flag dependencies are satisfied (e.g., `compute_cleavage` requires `relax_slabs`).

### Function Signature

```python
def resolve_workflow_preset(
    workflow_preset: Optional[str],
    relax_slabs: Optional[bool],
    compute_thermodynamics: Optional[bool],
    compute_cleavage: Optional[bool],
    compute_relaxation_energy: Optional[bool],
    compute_electronic_properties_bulk: Optional[bool],
    compute_electronic_properties_slabs: Optional[bool],
    run_aimd: Optional[bool],
    run_adsorption_energy: Optional[bool],
    # Validation inputs
    metal_name: Optional[str] = None,
    oxygen_name: Optional[str] = None,
    nonmetal_name: Optional[str] = None,
    miller_indices: Optional[list] = None,
    input_slabs: Optional[dict] = None,
    bands_parameters: Optional[dict] = None,
    bands_options: Optional[dict] = None,
    band_settings: Optional[dict] = None,
    aimd_sequence: Optional[list] = None,
    aimd_parameters: Optional[dict] = None,
    aimd_options: Optional[dict] = None,
    aimd_potential_mapping: Optional[dict] = None,
    aimd_kpoints_spacing: Optional[float] = None,
    slab_bands_parameters: Optional[dict] = None,
    slab_bands_options: Optional[dict] = None,
    slab_band_settings: Optional[dict] = None,
    slab_electronic_properties: Optional[dict] = None,
    adsorption_structures: Optional[dict] = None,
    adsorption_formulas: Optional[dict] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[str, dict]:
```

### Returns

- **Tuple of (resolved_preset_name, resolved_flags_dict)**:
  - `resolved_preset_name`: Name of the resolved preset (string)
  - `resolved_flags_dict`: Dictionary with all resolved flags:
    - `relax_slabs`: bool
    - `compute_thermodynamics`: bool
    - `compute_cleavage`: bool
    - `compute_relaxation_energy`: bool
    - `compute_electronic_properties_bulk`: bool
    - `compute_electronic_properties_slabs`: bool
    - `run_aimd`: bool
    - `run_adsorption_energy`: bool

### Raises

- **ValueError**: If:
  - `workflow_preset` is not a valid preset name
  - Required parameters for the preset are missing
  - Flag dependencies are not satisfied (errors, not warnings)

### Implementation Details

The function delegates to four helper functions from `teros.core.workflow_presets`:

1. `check_old_style_api()` - Emits deprecation warning if using old-style API
2. `resolve_preset()` - Loads preset and applies overrides
3. `validate_preset_inputs()` - Validates required parameters
4. `validate_flag_dependencies()` - Checks flag dependencies

### Usage Example

```python
from teros.core.stages import resolve_workflow_preset

# Basic usage with default preset
preset_name, flags = resolve_workflow_preset(
    workflow_preset=None,  # Uses default 'surface_thermodynamics'
    relax_slabs=None,
    compute_thermodynamics=None,
    compute_cleavage=None,
    compute_relaxation_energy=None,
    compute_electronic_properties_bulk=None,
    compute_electronic_properties_slabs=None,
    run_aimd=None,
    run_adsorption_energy=None,
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
)
# preset_name == 'surface_thermodynamics'
# flags['relax_slabs'] == True

# With explicit preset and override
preset_name, flags = resolve_workflow_preset(
    workflow_preset='surface_thermodynamics',
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=True,  # Override default (False)
    ...
)
# flags['compute_cleavage'] == True
```

### Files Created/Modified

- **Created**: `/home/thiagotd/git/PS-TEROS/teros/core/stages/preset_resolution.py` (~300 lines)
- **Modified**: `/home/thiagotd/git/PS-TEROS/teros/core/stages/__init__.py` (replaced placeholder with import)

### Verification

```bash
$ python3 -c "from teros.core.stages import resolve_workflow_preset; print('OK')"
OK
```

The function was tested with:
1. Default preset resolution (no preset specified)
2. Explicit preset with flag override
3. Both tests passed successfully

### Next Steps

The `resolve_workflow_preset()` function is now ready to be used in the refactored `build_core_workgraph()` function:

```python
# In build_core_workgraph():
from teros.core.stages import resolve_workflow_preset

resolved_preset_name, resolved_flags = resolve_workflow_preset(
    workflow_preset,
    relax_slabs,
    compute_thermodynamics,
    ...
)

# Extract resolved flags
relax_slabs = resolved_flags['relax_slabs']
compute_thermodynamics = resolved_flags['compute_thermodynamics']
# ... etc.
```
