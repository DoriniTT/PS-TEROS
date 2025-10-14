# PS-TEROS Workflow System Guide

## Overview

PS-TEROS provides two ways to configure workflows:

1. **Predefined Workflows** (`workflow_preset`): High-level presets for common use cases
2. **Custom Workflows**: Fine-grained control by setting individual flags

## Predefined Workflows

Use the `workflow_preset` parameter to quickly set up common workflows:

### Available Presets

#### 1. `surface_thermodynamics` (default)
Complete surface thermodynamics with relaxation.
- Relaxes slabs
- Computes surface energies as function of chemical potential γ(μ_O)
- **Cleavage and relaxation energies are optional** (set flags to enable)

**Required**: `metal_name`, `oxygen_name`

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    structures_dir='structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    # ... other parameters
)
```

#### 2. `surface_thermodynamics_unrelaxed`
Quick surface energy screening without slab relaxation.
- Uses unrelaxed slabs
- Fast initial assessment

**Required**: `metal_name`, `oxygen_name`

#### 3. `cleavage_only`
Calculate cleavage energies for slab pairs.
- Relaxes slabs
- Computes cleavage energies

#### 4. `relaxation_energy_only`
Calculate relaxation energies (E_relaxed - E_unrelaxed).
- Relaxes slabs
- Computes relaxation energies

#### 5. `bulk_only`
Bulk structure optimization only.
- No surface calculations

#### 6. `formation_enthalpy_only`
Formation enthalpy without surfaces.
- Relaxes bulk and references
- Computes formation enthalpy

**Required**: `metal_name`, `oxygen_name`

#### 7. `electronic_structure_bulk_only`
Electronic properties (DOS and bands) for bulk.
- Band structure calculation
- Density of states

**Required**: `bands_parameters`, `band_settings`

#### 8. `electronic_structure_slabs_only`
Electronic properties (DOS and bands) for slabs only.
- Relaxes slabs
- Computes band structure and DOS for slabs

**Required**: `slab_bands_parameters`, `slab_band_settings`

#### 9. `electronic_structure_bulk_and_slabs`
Electronic properties for both bulk and slabs.
- Complete electronic structure analysis

**Required**: `bands_parameters`, `band_settings`, `slab_bands_parameters`, `slab_band_settings`

#### 10. `aimd_only`
AIMD simulation on slabs.
- **Does NOT relax slabs before AIMD**
- Runs AIMD directly on generated slabs

**Required**: `aimd_sequence`, `aimd_parameters`

#### 11. `comprehensive`
Complete analysis with all features enabled.
- Thermodynamics + electronic properties + AIMD
- Most expensive workflow

## Custom Workflows

### Method 1: Set Individual Flags (No Preset)

For maximum flexibility, skip the `workflow_preset` and set flags individually:

```python
wg = build_core_workgraph(
    # NO workflow_preset parameter
    
    # Set individual flags
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

**See**: `step_10_custom_workflow.py`

### Method 2: Preset with Overrides

Start with a preset and override specific components:

```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    
    # Override preset defaults
    compute_cleavage=True,  # Add cleavage energies
    compute_relaxation_energy=True,  # Add relaxation energies
    
    # ... rest of parameters
)
```

**See**: `step_11_preset_with_overrides.py`

## Workflow Flags

All workflow flags (can be set individually or via preset):

- `relax_slabs`: Relax generated slabs
- `compute_thermodynamics`: Compute surface energies γ(μ)
- `compute_cleavage`: Compute cleavage energies
- `compute_relaxation_energy`: Compute relaxation energies
- `compute_electronic_properties_bulk`: Compute DOS/bands for bulk
- `compute_electronic_properties_slabs`: Compute DOS/bands for slabs
- `run_aimd`: Run AIMD simulations

## Examples by Use Case

### Surface Energy Calculations
- **Quick screening**: `surface_thermodynamics_unrelaxed`
- **Full analysis**: `surface_thermodynamics` + overrides

### Electronic Structure
- **Bulk only**: `electronic_structure_bulk_only`
- **Slabs only**: `electronic_structure_slabs_only`
- **Both**: `electronic_structure_bulk_and_slabs`

### Molecular Dynamics
- **AIMD**: `aimd_only`

### Cleavage/Relaxation
- **Cleavage**: `cleavage_only`
- **Relaxation**: `relaxation_energy_only`

### Complete Study
- **Everything**: `comprehensive`

## Step-by-Step Examples

| Step | File | Preset | Description |
|------|------|--------|-------------|
| 1 | `step_01_bulk_only.py` | `bulk_only` | Bulk relaxation |
| 2 | `step_02_formation_enthalpy.py` | `formation_enthalpy_only` | Formation enthalpy |
| 3 | `step_03_slabs_relaxation.py` | Custom | Slab relaxation |
| 4 | `step_04_cleavage_energy.py` | `cleavage_only` | Cleavage energies |
| 5 | `step_05_surface_thermodynamics.py` | `surface_thermodynamics` | Surface energies |
| 6 | `step_06_electronic_properties.py` | `electronic_structure_bulk_only` | Bulk DOS/bands |
| 7 | `step_07_aimd_simulation.py` | `aimd_only` | AIMD simulation |
| 8 | `step_08_electronic_structure_slabs.py` | `electronic_structure_slabs_only` | Slab DOS/bands |
| 9 | `step_09_electronic_structure_bulk_and_slabs.py` | `electronic_structure_bulk_and_slabs` | Full electronic structure |
| 10 | `step_10_custom_workflow.py` | None (custom) | Custom flag configuration |
| 11 | `step_11_preset_with_overrides.py` | `surface_thermodynamics` + overrides | Preset with additions |

## Key Changes from Previous Version

1. **Surface thermodynamics preset**: Cleavage and relaxation energies are now **optional** (not included by default)
2. **AIMD preset**: No longer relaxes slabs before AIMD (runs directly on generated slabs)
3. **New presets**: Added `electronic_structure_slabs_only` and `electronic_structure_bulk_and_slabs`
4. **More flexibility**: Easier to combine presets with custom overrides

## Listing Available Presets

```python
from teros.core.workflow_presets import list_workflow_presets

# Print all available presets
list_workflow_presets()

# Get detailed info about a specific preset
from teros.core.workflow_presets import get_preset_summary
print(get_preset_summary('surface_thermodynamics'))
```

## Tips

1. **Start with a preset** for common workflows
2. **Override specific flags** when you need custom behavior
3. **Use custom workflows** when combining unusual features
4. **Check validation** - presets will warn about missing required parameters
5. **Read the examples** - each step demonstrates a different pattern
