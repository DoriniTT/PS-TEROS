# Understanding the PS-TEROS Workflow System

## Overview

PS-TEROS uses a powerful three-tier workflow configuration system that balances simplicity with flexibility. This document explains how the system works and when to use each tier.

---

## The Three-Tier Architecture

### Tier 1: Named Workflow Presets (High-Level Convenience)

**What it is:** A single parameter that activates an entire workflow configuration.

**Example:**
```python
workflow_preset='surface_thermodynamics'
```

This one line automatically configures:
- ✓ Bulk relaxation
- ✓ Reference relaxations (metal, oxygen, optional nonmetal)
- ✓ Formation enthalpy calculation
- ✓ Slab generation and relaxation
- ✓ Surface energy calculations

**When to use:** For standard workflows and common use cases. This is the **recommended** approach for most users.

**Benefits:**
- Simple and clean
- Tested and validated
- Self-documenting code
- Reduces configuration errors

---

### Tier 2: Individual Component Flags (Fine-Grained Control)

**What it is:** Boolean flags that override preset defaults.

**Example:**
```python
workflow_preset='surface_thermodynamics',
compute_cleavage=True,              # Add cleavage energies
compute_relaxation_energy=False,    # Skip relaxation energies
```

**When to use:** When you need a preset workflow with minor modifications.

**Benefits:**
- Customize presets without starting from scratch
- Override specific features
- Maintain most preset benefits

**Available flags:**
- `relax_slabs`: Enable/disable slab relaxation
- `compute_thermodynamics`: Enable/disable surface energy calculations
- `compute_cleavage`: Enable/disable cleavage energy calculations
- `compute_relaxation_energy`: Enable/disable relaxation energy calculations
- `compute_electronic_properties_bulk`: Enable/disable bulk DOS/bands
- `compute_electronic_properties_slabs`: Enable/disable slab DOS/bands
- `run_aimd`: Enable/disable AIMD simulations

---

### Tier 3: Automatic Dependency Resolution (Smart Defaults)

**What it is:** The system automatically validates your configuration and warns about conflicts.

**Example scenario:**
```python
workflow_preset='surface_thermodynamics',
relax_slabs=False,  # You override this
compute_cleavage=True,  # But cleavage needs relaxed slabs!
```

**System response:**
```
⚠️ WARNING: compute_cleavage=True requires relax_slabs=True.
Cleavage energies will not be computed.
```

**Benefits:**
- Catches configuration errors early
- Prevents wasted computation time
- Clear error messages
- Automatic dependency checking

---

## How Preset Resolution Works

When you submit a workflow, the system follows this process:

```
1. Load preset configuration
   ↓
2. Apply user overrides (Tier 2 flags)
   ↓
3. Validate dependencies (Tier 3)
   ↓
4. Create workflow graph
```

### Visual Example

**Input:**
```python
workflow_preset='surface_thermodynamics',
compute_cleavage=False,  # User override
```

**Resolution:**
```
Preset 'surface_thermodynamics' defaults:
  relax_slabs = True
  compute_thermodynamics = True
  compute_cleavage = False  (OPTIONAL, default False)
  compute_relaxation_energy = False  (OPTIONAL, default False)
  compute_electronic_properties_bulk = False
  compute_electronic_properties_slabs = False
  run_aimd = False

User overrides:
  compute_cleavage = False  (no change from default)

Final configuration:
  relax_slabs = True
  compute_thermodynamics = True
  compute_cleavage = False
  compute_relaxation_energy = False
  compute_electronic_properties_bulk = False
  compute_electronic_properties_slabs = False
  run_aimd = False
```

---

## Decision Guide: Which Tier Should I Use?

### Use Tier 1 (Presets Only) When:
- ✓ Your use case matches a standard workflow
- ✓ You want the simplest configuration
- ✓ You're new to PS-TEROS
- ✓ You want tested, validated configurations

**Example:**
```python
workflow_preset='surface_thermodynamics',
# Just provide required parameters
```

---

### Use Tier 1 + Tier 2 (Preset with Overrides) When:
- ✓ You need a standard workflow with minor modifications
- ✓ You want to add optional features to a preset
- ✓ You want to skip certain calculations

**Example:**
```python
workflow_preset='surface_thermodynamics',
compute_cleavage=True,  # Add this optional feature
compute_relaxation_energy=True,  # Add this optional feature
```

---

### Use Tier 2 Only (Individual Flags) When:
- ✓ Your workflow doesn't match any preset
- ✓ You need a unique combination of features
- ✓ You're developing/testing new workflows

**Note:** You'll see a deprecation warning. Consider requesting a new preset if your use case is common.

**Example:**
```python
# No preset specified
relax_slabs=True,
compute_thermodynamics=False,
compute_cleavage=True,
# ... custom combination
```

---

## Understanding Preset Defaults

### Important: Optional vs. Required Components

Some presets have **optional** components that are disabled by default:

**Example: `surface_thermodynamics` preset**

**Always computed:**
- Bulk relaxation ✓
- Reference relaxations ✓
- Formation enthalpy ✓
- Slab generation and relaxation ✓
- Surface energies ✓

**Optional (disabled by default):**
- Cleavage energies ○ (add `compute_cleavage=True` to enable)
- Relaxation energies ○ (add `compute_relaxation_energy=True` to enable)

**Why optional?**
- Cleavage energies: Not always needed for thermodynamic analysis
- Relaxation energies: Additional computational cost, specific use case
- Flexibility: Users choose what they need

**How to enable:**
```python
workflow_preset='surface_thermodynamics',
compute_cleavage=True,              # Now it will compute cleavage
compute_relaxation_energy=True,     # Now it will compute relaxation
```

---

## Validation System

### What Gets Validated

**1. Required Parameters:**
```python
workflow_preset='surface_thermodynamics',
metal_name='Ag.cif',
# Missing oxygen_name!
```
**Error:**
```
ValueError: Preset 'surface_thermodynamics' requires parameter 'oxygen_name'
```

---

**2. Flag Dependencies:**
```python
workflow_preset='cleavage_only',
relax_slabs=False,  # This breaks cleavage calculation!
```
**Warning:**
```
⚠️ WARNING: compute_cleavage=True requires relax_slabs=True.
Cleavage energies will not be computed.
```

---

**3. Parameter Requirements:**
```python
workflow_preset='electronic_structure_bulk_only',
# Missing bands_parameters!
```
**Warning:**
```
⚠️ WARNING: compute_electronic_properties_bulk=True requires bands_parameters.
Electronic properties will not be computed.
```

---

## Common Patterns

### Pattern 1: Standard Workflow
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # ... required parameters
)
```
**Result:** Core thermodynamic calculations only (no optional features).

---

### Pattern 2: Add Optional Features
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,              # Enable optional
    compute_relaxation_energy=True,      # Enable optional
    # ... required parameters
)
```
**Result:** Core thermodynamics + cleavage + relaxation energies.

---

### Pattern 3: Remove Features
```python
wg = build_core_workgraph(
    workflow_preset='comprehensive',
    run_aimd=False,  # Skip AIMD
    # ... required parameters
)
```
**Result:** Everything except AIMD.

---

### Pattern 4: Inspect Before Running
```python
from teros.core.workflow_presets import resolve_preset, get_preset_summary

# See what a preset will do
print(get_preset_summary('surface_thermodynamics'))

# See final configuration with your overrides
preset_name, flags = resolve_preset(
    'surface_thermodynamics',
    compute_cleavage=True
)
print(f"Final flags: {flags}")
```
**Result:** Know exactly what will run before submitting.

---

## Troubleshooting

### Q: How do I know what a preset will compute?

**A:** Use `get_preset_summary()`:
```python
from teros.core import get_preset_summary
print(get_preset_summary('surface_thermodynamics'))
```

---

### Q: Why didn't my workflow compute what I expected?

**A:** Check the preset defaults. Some features are optional:
```python
from teros.core.workflow_presets import WORKFLOW_PRESETS
print(WORKFLOW_PRESETS['surface_thermodynamics']['flags'])
```

---

### Q: How do I see the final configuration before running?

**A:** Use `resolve_preset()`:
```python
from teros.core.workflow_presets import resolve_preset

preset_name, flags = resolve_preset(
    'surface_thermodynamics',
    compute_cleavage=True  # Your override
)
print(flags)
```

---

### Q: I'm getting a deprecation warning about using flags without presets

**A:** You're using Tier 2 without Tier 1. Either:
1. Switch to a preset: `workflow_preset='...'`
2. Request a new preset for your use case
3. Ignore the warning (not recommended)

---

### Q: Can I combine multiple presets?

**A:** No, only one preset can be active. However, you can:
1. Start with the closest preset
2. Override individual flags to customize it

---

## Best Practices

1. **Always start with a preset** - Only use individual flags for truly custom workflows
2. **Use `get_preset_summary()` before running** - Understand what will be computed
3. **Read validation warnings** - They save computation time
4. **Override minimally** - The fewer overrides, the clearer your intent
5. **Document custom configurations** - If using many overrides, explain why in comments

---

## Summary

The three-tier system gives you:
- **Simplicity** - One parameter for common workflows (Tier 1)
- **Flexibility** - Override any default (Tier 2)
- **Safety** - Automatic validation (Tier 3)

**Recommended workflow:**
1. Check available presets: `list_workflow_presets()`
2. Read preset details: `get_preset_summary('preset_name')`
3. Use preset with minimal overrides
4. Validate before running: `resolve_preset(...)`

---

## Related Documentation

- [Workflow Presets Guide](WORKFLOW_PRESETS_GUIDE.md) - Complete preset reference
- [Workflow Presets Examples](WORKFLOW_PRESETS_EXAMPLES.md) - Runnable code examples
- [Migration Guide](WORKFLOW_MIGRATION_GUIDE.md) - Updating from old API
