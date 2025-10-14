# PS-TEROS Workflow Quick Reference

## What Changed?

### 1. Surface Thermodynamics Preset
**Cleavage and relaxation energies are now OPTIONAL**

```python
# Default behavior (no cleavage or relaxation energies)
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # ... parameters
)

# To add cleavage and relaxation energies
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,  # ADD THIS
    compute_relaxation_energy=True,  # ADD THIS
    # ... parameters
)
```

### 2. AIMD Preset
**No longer relaxes slabs before AIMD**

```python
# Default behavior (no slab relaxation before AIMD)
wg = build_core_workgraph(
    workflow_preset='aimd_only',
    # ... parameters
)

# To relax slabs before AIMD (if needed)
wg = build_core_workgraph(
    workflow_preset='aimd_only',
    relax_slabs=True,  # ADD THIS
    # ... parameters
)
```

### 3. New Electronic Structure Presets

```python
# For slabs only
wg = build_core_workgraph(
    workflow_preset='electronic_structure_slabs_only',
    slab_bands_parameters=...,
    slab_band_settings=...,
    # ... parameters
)

# For both bulk and slabs
wg = build_core_workgraph(
    workflow_preset='electronic_structure_bulk_and_slabs',
    bands_parameters=...,
    band_settings=...,
    slab_bands_parameters=...,
    slab_band_settings=...,
    # ... parameters
)
```

## All Workflow Presets (11 total)

| Preset | What It Does | Use When |
|--------|-------------|----------|
| `surface_thermodynamics` | Surface energies γ(μ) | Main workflow |
| `surface_thermodynamics_unrelaxed` | Quick surface screening | Initial testing |
| `cleavage_only` | Cleavage energies | Material brittleness |
| `relaxation_energy_only` | Relaxation energies | Surface reconstruction |
| `bulk_only` | Bulk optimization | Testing setup |
| `formation_enthalpy_only` | Formation enthalpy | Stability analysis |
| `electronic_structure_bulk_only` | DOS/bands for bulk | Electronic properties |
| `electronic_structure_slabs_only` | DOS/bands for slabs | Surface electronics |
| `electronic_structure_bulk_and_slabs` | DOS/bands for both | Complete analysis |
| `aimd_only` | AIMD simulation | Dynamics studies |
| `comprehensive` | Everything enabled | Complete study |

## Two Ways to Configure Workflows

### Method 1: Use a Preset (Recommended)

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # ... rest of parameters
)
```

### Method 2: Set Individual Flags (Custom)

```python
wg = build_core_workgraph(
    # NO workflow_preset
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=True,
    compute_relaxation_energy=False,
    compute_electronic_properties_bulk=False,
    compute_electronic_properties_slabs=False,
    run_aimd=False,
    # ... rest of parameters
)
```

### Method 3: Preset + Overrides (Best of Both)

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,  # Override preset default
    # ... rest of parameters
)
```

## Examples

See `examples/step_by_step/`:

- **step_01-07**: Original examples (updated)
- **step_08**: Electronic structure for slabs only (NEW)
- **step_09**: Electronic structure for bulk and slabs (NEW)
- **step_10**: Custom workflow (no preset) (NEW)
- **step_11**: Preset with overrides (NEW)

## Full Documentation

- **Complete guide**: `examples/step_by_step/README_WORKFLOWS.md`
- **Step-by-step examples**: `examples/step_by_step/README.md`
- **Update summary**: `WORKFLOW_UPDATES_SUMMARY.md`

## List Available Presets

```python
from teros.core.workflow_presets import list_workflow_presets
list_workflow_presets()
```

## Get Preset Details

```python
from teros.core.workflow_presets import get_preset_summary
print(get_preset_summary('surface_thermodynamics'))
```
