# PS-TEROS Workflow Preset System - User Guide

## Overview

The PS-TEROS workflow preset system provides a simplified, three-tier approach to configuring complex computational workflows:

1. **Named Workflow Presets** (High-level convenience) - One parameter activates entire workflows
2. **Independent Component Flags** (Fine-grained control) - Override preset defaults as needed
3. **Automatic Dependency Resolution** (Smart defaults) - System validates and auto-enables dependencies

This guide explains how to use the preset system effectively.

---

## Quick Start

### Using a Preset (Recommended)

Instead of manually setting 7+ boolean flags, simply specify a preset:

```python
from teros.core.workgraph import build_core_workgraph

wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # One line activates full workflow!
    structures_dir='structures',
    bulk_name='ag3po4.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    bulk_potential_mapping={'Ag': 'Ag', 'P': 'P', 'O': 'O'},
    # ... other required parameters
)
```

### Listing Available Presets

To see all available presets:

```python
from teros.core import list_workflow_presets

list_workflow_presets()
```

Output shows each preset with description, use cases, and required parameters.

### Getting Preset Details

To see detailed configuration for a specific preset:

```python
from teros.core import get_preset_summary

print(get_preset_summary('surface_thermodynamics'))
```

---

## Available Workflow Presets

### 1. `surface_thermodynamics` (Default)

Complete surface thermodynamics workflow with relaxed slabs.

**Activates:**
- ✓ Bulk relaxation
- ✓ Reference relaxations (metal, oxygen, optional nonmetal)
- ✓ Formation enthalpy calculation
- ✓ Slab generation and relaxation
- ✓ Surface energy calculations with chemical potential sampling
- ○ Cleavage energy calculations (OPTIONAL - disabled by default, add `compute_cleavage=True` to enable)
- ○ Relaxation energy calculations (OPTIONAL - disabled by default, add `compute_relaxation_energy=True` to enable)

**Requirements:**
- `metal_name` (e.g., 'Ag.cif')
- `oxygen_name` (e.g., 'O2.cif')
- `miller_indices` or `input_slabs` (for slab generation)

**Use Cases:**
- Surface energy calculations
- Pourbaix-like stability analysis
- Complete thermodynamic characterization

**Note:** Cleavage and relaxation energies are optional features that can be enabled by setting their respective flags to `True`.

---

### 2. `surface_thermodynamics_unrelaxed`

Quick surface energy screening with unrelaxed slabs.

**Activates:**
- ✓ Bulk relaxation
- ✓ Reference relaxations
- ✓ Formation enthalpy calculation
- ✓ Slab generation (SCF only, no relaxation)
- ✓ Surface energy calculations

**Deactivates:**
- ✗ Slab relaxation
- ✗ Cleavage energies
- ✗ Relaxation energies

**Use Cases:**
- Quick screening of many terminations
- Testing slab generation parameters
- Low-cost initial assessment

---

### 3. `cleavage_only`

Cleavage energy calculations for complementary slab pairs.

**Activates:**
- ✓ Bulk relaxation
- ✓ Slab generation and relaxation
- ✓ Cleavage energy calculations

**Deactivates:**
- ✗ Reference relaxations
- ✗ Formation enthalpy
- ✗ Surface thermodynamics

**Use Cases:**
- Cleavage energy analysis
- Material brittleness assessment
- Surface creation energy studies

---

### 4. `relaxation_energy_only`

Calculate energy difference between unrelaxed and relaxed slabs.

**Activates:**
- ✓ Bulk relaxation
- ✓ Slab generation and relaxation
- ✓ Relaxation energy calculations

**Use Cases:**
- Surface reconstruction analysis
- Relaxation energy comparison
- Convergence studies

---

### 5. `bulk_only`

Bulk structure optimization only (no surfaces).

**Activates:**
- ✓ Bulk relaxation only

**Deactivates:**
- ✗ All surface-related calculations

**Use Cases:**
- Bulk structure optimization
- Testing calculation parameters
- Initial structure validation

---

### 6. `formation_enthalpy_only`

Formation enthalpy calculation without surface analysis.

**Activates:**
- ✓ Bulk relaxation
- ✓ Reference relaxations
- ✓ Formation enthalpy calculation

**Deactivates:**
- ✗ All surface-related calculations

**Requirements:**
- `metal_name`
- `oxygen_name`
- Optional: `nonmetal_name` (for ternary oxides)

**Use Cases:**
- Stability analysis
- Phase diagram construction
- Thermodynamic validation

---

### 7. `electronic_structure_bulk_only`

Electronic properties (DOS/bands) for bulk structure.

**Activates:**
- ✓ Bulk relaxation
- ✓ Electronic structure calculation (DOS and bands)

**Requirements:**
- `bands_parameters` - Dict with 'scf', 'bands', 'dos' INCAR parameters

**Use Cases:**
- Band structure analysis
- Density of states calculation
- Electronic property characterization

---

### 8. `electronic_structure_slabs_only`

Electronic properties (DOS/bands) for slabs only.

**Activates:**
- ✓ Bulk relaxation
- ✓ Slab generation and relaxation
- ✓ Electronic structure calculation for slabs (DOS and bands)

**Deactivates:**
- ✗ Bulk electronic properties
- ✗ Thermodynamics
- ✗ Cleavage/relaxation energies

**Requirements:**
- `slab_bands_parameters` - Dict with 'scf', 'bands', 'dos' INCAR parameters for slabs
- `slab_band_settings` - Band workflow settings for slabs
- `miller_indices` or `input_slabs`

**Use Cases:**
- Slab band structure calculation
- Slab density of states (DOS)
- Surface electronic property analysis

---

### 9. `electronic_structure_bulk_and_slabs`

Electronic properties (DOS/bands) for both bulk and slabs.

**Activates:**
- ✓ Bulk relaxation
- ✓ Slab generation and relaxation
- ✓ Electronic structure calculation for bulk (DOS and bands)
- ✓ Electronic structure calculation for slabs (DOS and bands)

**Deactivates:**
- ✗ Thermodynamics
- ✗ Cleavage/relaxation energies

**Requirements:**
- `bands_parameters` - Dict with 'scf', 'bands', 'dos' INCAR parameters for bulk
- `band_settings` - Band workflow settings for bulk
- `slab_bands_parameters` - Dict with 'scf', 'bands', 'dos' INCAR parameters for slabs
- `slab_band_settings` - Band workflow settings for slabs
- `miller_indices` or `input_slabs`

**Use Cases:**
- Complete electronic structure analysis
- Bulk vs surface electronic comparison
- Comprehensive band structure study

---

### 10. `aimd_only`

AIMD (Ab Initio Molecular Dynamics) simulation on slabs.

**Activates:**
- ✓ Bulk relaxation
- ✓ Slab generation (NO relaxation by default)
- ✓ AIMD simulation

**Deactivates:**
- ✗ Slab relaxation (add `relax_slabs=True` if you need slabs relaxed before AIMD)

**Requirements:**
- `aimd_sequence` - List of AIMD stages: `[{'temperature': K, 'steps': N}, ...]`
- `aimd_parameters` - AIMD INCAR parameters

**Use Cases:**
- Molecular dynamics simulation
- Temperature effects
- Dynamic property analysis

**Note:** By default, AIMD runs on unrelaxed slabs (generated from bulk). If you need relaxed slabs before AIMD, set `relax_slabs=True`.

---

### 11. `comprehensive`

Complete analysis with all features enabled.

**Activates:**
- ✓ Everything: thermodynamics, electronic properties, AIMD

**Requirements:**
- All parameters for thermodynamics, electronic properties, and AIMD

**Use Cases:**
- Complete material characterization
- Publication-ready comprehensive analysis
- Full property investigation

---

## Overriding Preset Defaults

You can override any preset default by explicitly setting individual flags:

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=False,  # Override: skip cleavage energies
    compute_relaxation_energy=False,  # Override: skip relaxation energies
    # ... rest of parameters
)
```

This gives you the convenience of presets with the flexibility of fine-grained control.

---

## Concurrency Control

All workflow builders accept a `max_concurrent_jobs` parameter to control how many VASP calculations run simultaneously:

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    max_concurrent_jobs=4,  # Default: max 4 VASP jobs at once
    # ... other parameters
)
```

**Values:**
- `1`: Serial mode (one VASP at a time)
- `4`: Default (moderate concurrency)
- `8+`: Higher concurrency for larger clusters
- `None`: Unlimited (full parallel)

See [CONCURRENCY_CONTROL.md](./CONCURRENCY_CONTROL.md) for details.

---

## Advanced Usage

### Custom Workflow Configuration

For unique workflows not covered by presets, you can still use individual flags:

```python
wg = build_core_workgraph(
    # No preset specified
    relax_slabs=True,
    compute_thermodynamics=False,
    compute_cleavage=True,
    compute_relaxation_energy=False,
    # ... rest of parameters
)
```

**Note:** A deprecation warning will be shown. Consider requesting a new preset if your use case is common.

### Checking Preset Configuration

Before running a workflow, you can inspect the resolved configuration:

```python
from teros.core.workflow_presets import resolve_preset

preset_name, flags = resolve_preset(
    'surface_thermodynamics',
    compute_cleavage=False  # Your override
)

print(f"Preset: {preset_name}")
print(f"Flags: {flags}")
```

---

## Parameter Requirements by Preset

| Preset | metal_name | oxygen_name | nonmetal_name | miller_indices | bands_params | slab_bands_params | aimd_params |
|--------|------------|-------------|---------------|----------------|--------------|-------------------|-------------|
| `surface_thermodynamics` | ✓ | ✓ | ○ | ○ | ✗ | ✗ | ✗ |
| `surface_thermodynamics_unrelaxed` | ✓ | ✓ | ○ | ○ | ✗ | ✗ | ✗ |
| `cleavage_only` | ✗ | ✗ | ✗ | ○ | ✗ | ✗ | ✗ |
| `relaxation_energy_only` | ✗ | ✗ | ✗ | ○ | ✗ | ✗ | ✗ |
| `bulk_only` | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ |
| `formation_enthalpy_only` | ✓ | ✓ | ○ | ✗ | ✗ | ✗ | ✗ |
| `electronic_structure_bulk_only` | ✗ | ✗ | ✗ | ✗ | ✓ | ✗ | ✗ |
| `electronic_structure_slabs_only` | ✗ | ✗ | ✗ | ○ | ✗ | ✓ | ✗ |
| `electronic_structure_bulk_and_slabs` | ✗ | ✗ | ✗ | ○ | ✓ | ✓ | ✗ |
| `aimd_only` | ✗ | ✗ | ✗ | ○ | ✗ | ✗ | ✓ |
| `comprehensive` | ✓ | ✓ | ○ | ○ | ✓ | ✓ | ✓ |

**Legend:**
- ✓ = Required
- ○ = Optional but recommended
- ✗ = Not needed

---

## Validation and Error Handling

The preset system validates your configuration and provides clear error messages:

### Missing Required Parameters

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    metal_name='Ag.cif',
    # Missing oxygen_name!
)
```

**Error:**
```
ValueError: Preset validation failed:
Preset 'surface_thermodynamics' requires parameter 'oxygen_name' but it was not provided or is None
```

### Flag Dependencies

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    relax_slabs=False,  # Override
    # This creates inconsistency!
)
```

**Warning:**
```
⚠️  WARNING: compute_cleavage=True requires relax_slabs=True.
Cleavage energies will not be computed.
```

The workflow will still run, but affected calculations will be skipped.

---

## Migration from Old API

### Old Style (Deprecated)

```python
wg = build_core_workgraph(
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=True,
    compute_relaxation_energy=True,
    compute_electronic_properties_bulk=False,
    compute_electronic_properties_slabs=False,
    run_aimd=False,
    # ... rest
)
```

### New Style (Recommended)

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # ... rest (all flags set automatically)
)
```

---

## Best Practices

1. **Always use presets** for common workflows - they're tested and validated
2. **Override sparingly** - only when you have a specific reason
3. **Check preset summaries** before using them: `get_preset_summary('preset_name')`
4. **Read validation warnings** - they help catch configuration errors early
5. **Request new presets** - if your use case isn't covered, suggest a new preset

---

## Troubleshooting

### Q: Which preset should I use?

**A:** Start with `surface_thermodynamics` for most surface calculations. Use `list_workflow_presets()` to explore options.

### Q: Can I combine presets?

**A:** No, only one preset can be active. However, you can start with a preset and override individual flags.

### Q: My workflow isn't covered by any preset. What do I do?

**A:** Use individual flags (you'll see a deprecation warning). Consider requesting a new preset for common use cases.

### Q: How do I disable a feature from a preset?

**A:** Override the specific flag:
```python
workflow_preset='surface_thermodynamics',
compute_cleavage=False,  # Disable this feature
```

### Q: What happens if I don't specify a preset?

**A:** The system defaults to `surface_thermodynamics`. If you set individual flags without a preset, you'll see a deprecation warning.

---

## Summary

The workflow preset system simplifies PS-TEROS usage while maintaining full flexibility:

- **Simple:** One parameter for common workflows
- **Flexible:** Override any default
- **Safe:** Automatic validation and dependency checking
- **Clear:** Helpful error messages and warnings

For examples, see [WORKFLOW_PRESETS_EXAMPLES.md](WORKFLOW_PRESETS_EXAMPLES.md).

For implementation details, see the [Implementation Guide](../WORKFLOW_PRESET_IMPLEMENTATION_GUIDE.md).
