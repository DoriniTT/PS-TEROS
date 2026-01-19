# Slab Electronic Properties Stage Binding

## Summary
Replaced inline slab electronic properties code in `build_core_workgraph()` with a call to `add_slab_electronic_properties_stage()`.

## Changes Made

### File: `teros/core/workgraph.py`

**Lines replaced:** 1564-1612 (original line numbers, now at ~1736-1766)

**Before:** 49 lines of inline code including:
- Import of `calculate_electronic_properties_slabs_scatter`
- Calculation of default parameters, options, and settings
- Code loading via `load_code()`
- Computation of slab parameters/options/mapping with fallbacks
- Conditional source selection based on restart mode
- Task creation with `wg.add_task()`
- Output socket connections
- Logger call

**After:** 22 lines with stage function call:
```python
# Stage 7: Add slab electronic properties calculation if requested
# Currently only supported for input_slabs mode (including restart)
if (
    compute_electronic_properties_slabs
    and relax_slabs
    and slab_electronic_properties
    and use_input_slabs
):
    # Get slab potential mapping with fallback
    slab_pot_map = (
        slab_potential_mapping
        if slab_potential_mapping is not None
        else bulk_potential_mapping
    )

    add_slab_electronic_properties_stage(
        wg=wg,
        code_label=code_label,
        potential_family=potential_family,
        bulk_options=bulk_options,
        slab_potential_mapping=slab_pot_map,
        slab_bands_parameters=slab_bands_parameters,
        slab_bands_options=slab_bands_options,
        slab_band_settings=slab_band_settings,
        slab_electronic_properties=slab_electronic_properties,
        max_concurrent_jobs=max_concurrent_jobs,
        use_restart_mode=restart_folders is not None,
        logger=logger,
        scatter_task=scatter_task,
        collector=collector,
    )
```

## Behavioral Equivalence

| Aspect | Before | After |
|--------|--------|-------|
| Conditional check | Same 4-part condition | Same 4-part condition |
| Potential mapping | Fallback to `bulk_potential_mapping` | Fallback computed before call |
| Restart mode detection | `restart_folders is not None` | Passed as `use_restart_mode` boolean |
| Task outputs | 4 outputs connected | Same 4 outputs (in stage function) |
| Logger message | `"✓ Slab electronic properties..."` | Same message (in stage function) |

## Import Required (deferred)
```python
from teros.core.stages.electronic_properties import add_slab_electronic_properties_stage
```

Note: Import will be added in a later consolidation step.

## Line Count Reduction
- Before: 49 lines of inline code
- After: 22 lines (including blank lines and comments)
- Net reduction: ~27 lines (55% reduction)
