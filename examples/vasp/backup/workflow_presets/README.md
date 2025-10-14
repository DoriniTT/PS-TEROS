# PS-TEROS Workflow Preset Examples

This directory contains example scripts demonstrating the workflow preset system.

## Quick Start

1. **List available presets:**
   ```bash
   python -c "from teros.core import list_workflow_presets; list_workflow_presets()"
   ```

2. **Run a simple example:**
   ```bash
   python example_bulk_only.py
   ```

## Available Examples

### `example_surface_thermo.py`
Complete surface thermodynamics workflow using the default `surface_thermodynamics` preset.

**Calculates:**
- Bulk and reference relaxations
- Formation enthalpy
- Slab generation and relaxation
- Surface energies
- Cleavage energies
- Relaxation energies

**Good for:** Complete surface characterization

---

### `example_bulk_only.py`
Simple bulk structure optimization.

**Calculates:**
- Bulk relaxation only

**Good for:** Testing parameters, initial validation

---

### `example_aimd.py`
AIMD (Ab Initio Molecular Dynamics) simulation on surface slabs.

**Calculates:**
- Bulk relaxation
- Slab generation and relaxation
- AIMD simulation (equilibration + production)

**Good for:** Dynamic properties, temperature effects

---

### `example_custom_workflow.py`
Demonstrates how to override preset defaults for custom workflows.

**Shows:**
- Starting with a preset
- Overriding specific components
- Interactive preset exploration

**Good for:** Learning how to customize workflows

---

## Before Running

1. **Update configuration** in each example:
   - Change `STRUCTURES_DIR` to your structures directory
   - Change `CODE_LABEL` to your VASP code label
   - Adjust VASP parameters and scheduler options as needed

2. **Ensure AiiDA is running:**
   ```bash
   verdi status
   verdi daemon restart
   ```

3. **Check available presets:**
   ```python
   from teros.core import list_workflow_presets
   list_workflow_presets()
   ```

## Workflow Presets Summary

| Preset | Description | Required Parameters |
|--------|-------------|-------------------|
| `surface_thermodynamics` | Complete surface thermodynamics | metal_name, oxygen_name |
| `surface_thermodynamics_unrelaxed` | Quick screening (unrelaxed) | metal_name, oxygen_name |
| `cleavage_only` | Cleavage energies only | — |
| `relaxation_energy_only` | Relaxation energies only | — |
| `bulk_only` | Bulk relaxation only | — |
| `formation_enthalpy_only` | Formation enthalpy only | metal_name, oxygen_name |
| `electronic_structure_bulk_only` | DOS/bands for bulk | bands_parameters |
| `aimd_only` | AIMD simulation | aimd_sequence, aimd_parameters |
| `comprehensive` | Everything enabled | All of the above |

## Documentation

- **User Guide:** `docs/WORKFLOW_PRESETS_GUIDE.md`
- **More Examples:** `docs/WORKFLOW_PRESETS_EXAMPLES.md`
- **Implementation:** `WORKFLOW_PRESET_IMPLEMENTATION_GUIDE.md`

## Tips

1. **Start simple:** Begin with `bulk_only` or `formation_enthalpy_only`
2. **Use presets:** They handle dependencies automatically
3. **Override when needed:** Customize specific components
4. **Check validation:** Preset system validates your configuration
5. **Monitor workflows:** Use `verdi process status <PK>`

## Getting Help

```python
# List all presets with descriptions
from teros.core import list_workflow_presets
list_workflow_presets()

# Get detailed info about a specific preset
from teros.core import get_preset_summary
print(get_preset_summary('surface_thermodynamics'))

# Get preset configuration
from teros.core import get_preset_config
config = get_preset_config('bulk_only')
print(config)
```
