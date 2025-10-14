# PS-TEROS Workflow Migration Guide

## Overview

This guide helps you update existing PS-TEROS scripts to work with the updated workflow preset system introduced in October 2025.

**Good news:** All changes are **backward compatible**. Your existing scripts will still work, though you may see deprecation warnings.

---

## What Changed?

### Summary of Changes

| Component | Previous Behavior | New Behavior | Action Required |
|-----------|------------------|--------------|-----------------|
| `surface_thermodynamics` preset | Always computed cleavage and relaxation energies | Cleavage and relaxation are **optional** (disabled by default) | Add explicit flags if you need these features |
| `aimd_only` preset | Automatically relaxed slabs before AIMD | Only generates slabs (no relaxation) | Add `relax_slabs=True` if you need relaxation first |
| New presets | N/A | Added `electronic_structure_slabs_only` and `electronic_structure_bulk_and_slabs` | Use new presets for electronic structure calculations |

---

## Breaking vs. Non-Breaking Changes

### ✅ Backward Compatible (Non-Breaking)

These changes **will not break** existing scripts:
- All presets still exist with the same names
- All parameters are still accepted
- Old flag-based API still works (with deprecation warning)
- Default preset is still `surface_thermodynamics`

### ⚠️ Behavioral Changes (May Affect Results)

These changes **may produce different outputs** than before:

**1. `surface_thermodynamics` preset:**
   - **Before:** Automatically computed cleavage and relaxation energies
   - **Now:** These are optional, disabled by default
   - **Impact:** Workflows using this preset will no longer compute cleavage/relaxation unless explicitly enabled

**2. `aimd_only` preset:**
   - **Before:** Relaxed slabs before running AIMD
   - **Now:** Only generates slabs (no relaxation step)
   - **Impact:** AIMD will run on unrelaxed slabs unless you add `relax_slabs=True`

---

## Do I Need to Update My Scripts?

### Decision Flowchart

```
Do you use workflow_preset='surface_thermodynamics'?
├─ YES → Do you need cleavage or relaxation energies?
│         ├─ YES → UPDATE REQUIRED (see Section A)
│         └─ NO → No action needed
│
└─ NO → Do you use workflow_preset='aimd_only'?
          ├─ YES → Do you need slabs relaxed before AIMD?
          │         ├─ YES → UPDATE REQUIRED (see Section B)
          │         └─ NO → No action needed
          │
          └─ NO → Do you use individual flags without a preset?
                    ├─ YES → Consider using presets (see Section C)
                    └─ NO → No action needed
```

---

## Migration Scenarios

### Section A: You Use `surface_thermodynamics` Preset

#### Scenario A1: You Need Cleavage Energies

**Old script (implicit behavior):**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    structures_dir='structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    # ... other parameters
)
# Previously, this automatically computed cleavage energies
```

**New script (explicit flags):**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,  # ← ADD THIS
    structures_dir='structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    # ... other parameters
)
```

---

#### Scenario A2: You Need Relaxation Energies

**Old script (implicit behavior):**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # ... parameters
)
# Previously, this automatically computed relaxation energies
```

**New script (explicit flags):**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_relaxation_energy=True,  # ← ADD THIS
    # ... parameters
)
```

---

#### Scenario A3: You Need Both

**New script:**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,              # ← ADD THIS
    compute_relaxation_energy=True,      # ← ADD THIS
    # ... parameters
)
```

---

#### Scenario A4: You Don't Need Either (No Action)

**Your script:**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # ... parameters
)
```

**Result:** ✅ No changes needed. Script works as-is and produces the same core results (surface energies, formation enthalpy, etc.)

---

### Section B: You Use `aimd_only` Preset

#### Scenario B1: You Need Relaxed Slabs for AIMD

**Old script (implicit behavior):**
```python
wg = build_core_workgraph(
    workflow_preset='aimd_only',
    # ... parameters
    aimd_sequence=[{'temperature': 300, 'steps': 1000}],
)
# Previously, this relaxed slabs before AIMD
```

**New script (explicit flag):**
```python
wg = build_core_workgraph(
    workflow_preset='aimd_only',
    relax_slabs=True,  # ← ADD THIS if you want relaxation
    # ... parameters
    aimd_sequence=[{'temperature': 300, 'steps': 1000}],
)
```

---

#### Scenario B2: You Want AIMD on Unrelaxed Slabs (No Action)

**Your script:**
```python
wg = build_core_workgraph(
    workflow_preset='aimd_only',
    # ... parameters
)
```

**Result:** ✅ No changes needed. This is the new default behavior.

---

### Section C: You Use Individual Flags Without Presets

#### Scenario C1: Migrate to Presets (Recommended)

**Old script:**
```python
wg = build_core_workgraph(
    # No preset
    relax_slabs=True,
    compute_thermodynamics=True,
    compute_cleavage=True,
    compute_relaxation_energy=False,
    compute_electronic_properties_bulk=False,
    compute_electronic_properties_slabs=False,
    run_aimd=False,
    # ... parameters
)
# Warning: Deprecation warning shown
```

**New script:**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # ← USE PRESET
    compute_cleavage=True,  # ← Only override what's different
    # ... parameters
)
# Cleaner, self-documenting, no warning
```

---

#### Scenario C2: Keep Using Flags (Not Recommended)

**Your script:**
```python
wg = build_core_workgraph(
    relax_slabs=True,
    compute_thermodynamics=False,
    compute_cleavage=True,
    # ... custom combination
)
# Still works, but shows deprecation warning
```

**Result:** Script still works. You'll see:
```
DeprecationWarning: Calling build_core_workgraph() with explicit boolean flags
without specifying workflow_preset is deprecated. Please use workflow_preset
parameter for better clarity and maintainability.
```

**Recommendation:** Consider requesting a new preset if your use case is common.

---

## Testing Your Migration

### Step 1: Check Current Behavior

Before updating, document what your current workflow produces:

```bash
# Run existing script
python my_workflow.py

# Check outputs
verdi process show <PK>

# Note which outputs exist:
# - cleavage_energies?
# - relaxation_energies?
# - slab_energies (relaxed)?
```

---

### Step 2: Update Script

Apply the appropriate changes from Sections A, B, or C above.

---

### Step 3: Verify Configuration

Before submitting, check what the new configuration will do:

```python
from teros.core.workflow_presets import resolve_preset

preset_name, flags = resolve_preset(
    'surface_thermodynamics',
    compute_cleavage=True,  # Your overrides
)

print(f"Preset: {preset_name}")
print(f"Flags: {flags}")
```

---

### Step 4: Test Run

Submit a test workflow and verify outputs match expectations:

```bash
python updated_workflow.py
verdi process show <NEW_PK>

# Verify same outputs as Step 1
```

---

## Common Migration Mistakes

### Mistake 1: Assuming Default Behavior Hasn't Changed

**Problem:**
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    # Expecting cleavage energies...
)
```

**Solution:** Explicitly enable optional features:
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',
    compute_cleavage=True,  # Make expectations explicit
)
```

---

### Mistake 2: Not Reading Preset Defaults

**Problem:** Using a preset without understanding what it does.

**Solution:** Always check before using:
```python
from teros.core import get_preset_summary
print(get_preset_summary('surface_thermodynamics'))
```

---

### Mistake 3: Ignoring Validation Warnings

**Problem:** Dismissing warnings during workflow submission.

**Solution:** Read warnings carefully - they indicate configuration issues:
```
⚠️ WARNING: compute_cleavage=True requires relax_slabs=True.
```

---

## Migration Checklist

Use this checklist when updating scripts:

- [ ] Identify which preset(s) you use
- [ ] Check if preset behavior changed (see "What Changed?" section)
- [ ] Document current workflow outputs
- [ ] Apply appropriate migration changes
- [ ] Verify new configuration with `resolve_preset()`
- [ ] Test updated script with small example
- [ ] Compare outputs before/after
- [ ] Update any automation/CI scripts
- [ ] Update documentation/comments in code

---

## Getting Help

### Q: How do I know which preset I'm using?

**A:** Check your script for the `workflow_preset` parameter:
```python
wg = build_core_workgraph(
    workflow_preset='surface_thermodynamics',  # ← Here
    # ...
)
```

If you don't see `workflow_preset`, you're using individual flags (Section C).

---

### Q: My workflow used to produce X, but now it doesn't

**A:** Check if X is an optional feature:
- Cleavage energies: Add `compute_cleavage=True`
- Relaxation energies: Add `compute_relaxation_energy=True`
- Relaxed slabs before AIMD: Add `relax_slabs=True` to `aimd_only` preset

---

### Q: Can I keep using the old API?

**A:** Yes, but not recommended. You'll see deprecation warnings. The preset system is clearer and more maintainable.

---

### Q: What if my use case doesn't match any preset?

**A:**
1. Use individual flags (you'll see a deprecation warning)
2. Request a new preset by opening a GitHub issue
3. Check if a combination of preset + overrides works

---

## New Features You Can Use

While migrating, consider taking advantage of new features:

### New Preset: `electronic_structure_slabs_only`

Calculate DOS and band structure for slabs:
```python
wg = build_core_workgraph(
    workflow_preset='electronic_structure_slabs_only',
    # ... parameters
    slab_bands_parameters=...,
)
```

---

### New Preset: `electronic_structure_bulk_and_slabs`

Calculate electronic properties for both bulk and slabs:
```python
wg = build_core_workgraph(
    workflow_preset='electronic_structure_bulk_and_slabs',
    # ... parameters
    bands_parameters=...,
    slab_bands_parameters=...,
)
```

---

## Summary

**Key takeaways:**
1. Changes are backward compatible
2. `surface_thermodynamics`: cleavage/relaxation now optional (add explicit flags if needed)
3. `aimd_only`: no automatic relaxation (add `relax_slabs=True` if needed)
4. Consider migrating from flags to presets
5. Always verify configuration before running expensive workflows

**Recommended workflow:**
1. Check current outputs
2. Apply migration changes
3. Verify with `resolve_preset()`
4. Test with small example
5. Update production scripts

---

## Related Documentation

- [Workflow System Explained](WORKFLOW_SYSTEM_EXPLAINED.md) - Understanding the three-tier system
- [Workflow Presets Guide](WORKFLOW_PRESETS_GUIDE.md) - Complete preset reference
- [Workflow Presets Examples](WORKFLOW_PRESETS_EXAMPLES.md) - Code examples
