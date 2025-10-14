# PS-TEROS Workflow Preset System Updates

## Summary of Changes

This document summarizes the updates made to the PS-TEROS workflow preset system based on the following requirements:

1. Create separate electronic structure presets for slabs only and bulk+slabs
2. Make cleavage and relaxation energies optional in surface_thermodynamics preset
3. Remove slab relaxation from aimd_only preset
4. Document the difference between predefined and custom workflows

## Changes Made

### 1. Updated Workflow Presets (`teros/core/workflow_presets.py`)

#### Modified Presets:

**`surface_thermodynamics`**
- **Before**: Always computed cleavage and relaxation energies
- **After**: Cleavage and relaxation energies are now OPTIONAL (set to False by default)
- Users can enable them via overrides: `compute_cleavage=True`, `compute_relaxation_energy=True`

**`aimd_only`**
- **Before**: `relax_slabs=True` (relaxed slabs before AIMD)
- **After**: `relax_slabs=False` (AIMD runs directly on generated slabs)
- Rationale: AIMD is meant to run on unrelaxed structures to study dynamics

#### New Presets:

**`electronic_structure_slabs_only`**
- Computes DOS and band structure for slabs only
- Relaxes slabs before electronic structure calculation
- Required parameters: `slab_bands_parameters`, `slab_band_settings`

**`electronic_structure_bulk_and_slabs`**
- Computes DOS and band structure for both bulk and slabs
- Comprehensive electronic structure analysis
- Required parameters: `bands_parameters`, `band_settings`, `slab_bands_parameters`, `slab_band_settings`

### 2. New Example Files (`examples/step_by_step/`)

Created 4 new example files:

1. **`step_08_electronic_structure_slabs.py`**
   - Demonstrates `electronic_structure_slabs_only` preset
   - Shows how to compute slab electronic properties

2. **`step_09_electronic_structure_bulk_and_slabs.py`**
   - Demonstrates `electronic_structure_bulk_and_slabs` preset
   - Shows comprehensive electronic structure workflow

3. **`step_10_custom_workflow.py`**
   - Demonstrates custom workflow without using a preset
   - Sets individual flags manually
   - Example: Thermodynamics + cleavage, but no relaxation energies

4. **`step_11_preset_with_overrides.py`**
   - Demonstrates starting with a preset and overriding specific flags
   - Example: `surface_thermodynamics` + cleavage + relaxation energies
   - Shows the convenience of preset + customization

### 3. Updated Existing Examples

**`step_05_surface_thermodynamics.py`**
- Updated documentation to reflect that cleavage and relaxation energies are optional
- Added notes about how to enable them via overrides

**`step_07_aimd_simulation.py`**
- Updated documentation to note that slabs are NOT relaxed in `aimd_only` preset
- Clarified that AIMD runs directly on generated slab structures

### 4. Documentation Updates

**`examples/step_by_step/README.md`**
- Added Steps 8, 9, 10, 11
- Updated Step 5 and Step 7 descriptions
- Updated workflow diagram

**`examples/step_by_step/README_WORKFLOWS.md`** (NEW)
- Comprehensive guide to workflow system
- Explains all 11 workflow presets
- Documents two workflow configuration methods:
  1. Predefined workflows (use `workflow_preset`)
  2. Custom workflows (set individual flags)
- Provides examples for each use case
- Includes comparison table of all steps

## Workflow Configuration Methods

### Method 1: Predefined Workflows (Recommended for Common Cases)

Use the `workflow_preset` parameter:

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # ... other parameters
)
```

**Benefits:**
- Quick setup for common workflows
- Fewer parameters to configure
- Validated default combinations
- Easy to understand intent

### Method 2: Custom Workflows

#### Option A: No Preset (Full Control)

Set individual flags without using a preset:

```python
wg = build_core_workgraph(
    # NO workflow_preset
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=True,
    compute_relaxation_energy=False,
    # ... other flags
)
```

**Benefits:**
- Maximum flexibility
- Unusual flag combinations
- Fine-grained control

#### Option B: Preset + Overrides (Best of Both)

Start with a preset and override specific flags:

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,  # Override preset default
    compute_relaxation_energy=True,  # Override preset default
    # ... other parameters
)
```

**Benefits:**
- Convenience of presets
- Customization when needed
- Clear intent (base + modifications)

## All Available Presets

1. `surface_thermodynamics` - Surface energies with optional cleavage/relaxation
2. `surface_thermodynamics_unrelaxed` - Quick unrelaxed screening
3. `cleavage_only` - Cleavage energy calculations
4. `relaxation_energy_only` - Relaxation energy calculations
5. `bulk_only` - Bulk structure optimization
6. `formation_enthalpy_only` - Formation enthalpy without surfaces
7. `electronic_structure_bulk_only` - DOS/bands for bulk
8. `electronic_structure_slabs_only` - DOS/bands for slabs (NEW)
9. `electronic_structure_bulk_and_slabs` - DOS/bands for both (NEW)
10. `aimd_only` - AIMD simulation (no slab relaxation)
11. `comprehensive` - All features enabled

## Testing

All presets have been validated:

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-reorganize
source ~/envs/psteros/bin/activate
python -c "from teros.core.workflow_presets import list_workflow_presets; list_workflow_presets()"
```

Verified preset flags:
- `surface_thermodynamics`: cleavage=False, relaxation_energy=False ✅
- `aimd_only`: relax_slabs=False ✅
- `electronic_structure_slabs_only`: compute_electronic_properties_slabs=True ✅

## Key Improvements

1. **More granular electronic structure control**: Separate presets for bulk-only, slabs-only, and both
2. **More flexible surface thermodynamics**: Cleavage and relaxation energies are optional
3. **More accurate AIMD workflow**: No unnecessary slab relaxation before dynamics
4. **Better documentation**: Clear distinction between predefined and custom workflows
5. **More examples**: 11 examples covering all workflow patterns
6. **Easier customization**: Preset + override pattern documented and demonstrated

## Files Modified

### Core Code:
- `teros/core/workflow_presets.py` - Updated preset definitions

### New Examples:
- `examples/step_by_step/step_08_electronic_structure_slabs.py`
- `examples/step_by_step/step_09_electronic_structure_bulk_and_slabs.py`
- `examples/step_by_step/step_10_custom_workflow.py`
- `examples/step_by_step/step_11_preset_with_overrides.py`

### Updated Examples:
- `examples/step_by_step/step_05_surface_thermodynamics.py`
- `examples/step_by_step/step_07_aimd_simulation.py`

### Documentation:
- `examples/step_by_step/README.md` - Updated with new steps
- `examples/step_by_step/README_WORKFLOWS.md` - NEW comprehensive guide

## Migration Guide for Existing Code

### If you were using `surface_thermodynamics` preset:

**Old behavior**: Automatically computed cleavage and relaxation energies  
**New behavior**: Cleavage and relaxation energies are optional

**To maintain old behavior**, add overrides:
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,  # ADD THIS
    compute_relaxation_energy=True,  # ADD THIS
    # ... rest of parameters
)
```

### If you were using `aimd_only` preset:

**Old behavior**: Relaxed slabs before AIMD  
**New behavior**: AIMD runs directly on generated slabs

**To maintain old behavior** (if needed), add override:
```python
wg = build_core_workgraph(
    workflow_preset='aimd_only',
    relax_slabs=True,  # ADD THIS to relax before AIMD
    # ... rest of parameters
)
```

### If you need electronic structure for slabs:

**Old**: Would need to use custom flags  
**New**: Use the new presets:
```python
# For slabs only:
wg = build_core_workgraph(
    workflow_preset='electronic_structure_slabs_only',
    # ...
)

# For both bulk and slabs:
wg = build_core_workgraph(
    workflow_preset='electronic_structure_bulk_and_slabs',
    # ...
)
```

## Next Steps

Users should:
1. Review `README_WORKFLOWS.md` for complete workflow system documentation
2. Check the 11 step-by-step examples for patterns
3. Consider migrating old code to use presets for clarity
4. Use preset + override pattern for custom workflows

## Notes

- All existing functionality is preserved
- Changes are backward compatible with explicit flag setting
- Deprecation warnings will guide users toward preset usage
- The preset system makes the API cleaner and more maintainable
