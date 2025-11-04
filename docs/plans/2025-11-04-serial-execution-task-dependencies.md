# Serial Execution with Task Dependencies

**Date:** 2025-11-04
**Status:** Approved
**Feature:** Hybrid serial execution control using task dependencies and max_concurrent_jobs

---

## Problem Statement

The current `max_concurrent_jobs` parameter controls global concurrency but does not enforce execution order between workflow stages. When `max_concurrent_jobs=1` on cluster03 (no queue system):

**Current behavior (undesired):**
```
Bulk [done] → Metal [waiting]
           ↓
           Slab_0 [running]  ← Problem: Slab launches before references finish
```

**Desired behavior:**
```
Bulk → Metal → Oxygen → SCF(term_0) → SCF(term_1) → Relax(term_0) → Relax(term_1)
     (strict serial execution, one VASP at a time globally)
```

## Requirements

1. **Fully serial execution** when `max_concurrent_jobs=1`
   - Only 1 VASP calculation running at ANY time across ENTIRE workflow

2. **Reference calculations in sequence:**
   - Bulk → Metal → Nonmetal (if ternary) → Oxygen

3. **Stage separation:**
   - All references must complete before any slab calculation starts

4. **Slab stage dependencies:**
   - When `compute_relaxation_energy=True`: SCF calculations complete before relaxations start
   - Chain: references → SCF → relaxation

5. **Apply to all PS-TEROS modules:**
   - Core workflows
   - AIMD module
   - Surface hydroxylation module
   - Custom calculation module

---

## Architecture

### Hybrid Approach

Combine two complementary mechanisms:

1. **Task dependency operators (`>>`):** Enforce execution order
2. **Global `max_number_jobs`:** Enforce concurrency limit

**Why both?**
- `>>` ensures references complete before slabs start
- `max_number_jobs` prevents multiple calculations in parallel
- Together: strict serial execution

### Dependency Chain Design

```python
# Reference calculations (serial)
bulk_vasp >> metal_vasp >> nonmetal_vasp >> oxygen_vasp

# Stage separation (references → slabs)
group(bulk_vasp, metal_vasp, oxygen_vasp) >> slab_stage

# Slab stage (SCF → relaxation)
scf_outputs >> relaxation_outputs
```

**Complete execution flow:**
```
Bulk → Metal → Nonmetal → Oxygen → group(all_references)
                                         ↓
                    SCF(term_0) → SCF(term_1) → ... → SCF(term_N)
                                                           ↓
                              Relax(term_0) → Relax(term_1) → ... → Relax(term_N)
```

---

## Implementation

### Module: `teros/core/workgraph.py`

**Changes required:**

1. **Import group utility:**
   ```python
   from aiida_workgraph import task, WorkGraph, dynamic, namespace, group
   ```

2. **Chain reference calculations:**
   ```python
   if compute_formation_enthalpy:
       # Create tasks
       metal_vasp = VaspTask(...)
       oxygen_vasp = VaspTask(...)

       # Chain dependencies
       bulk_vasp >> metal_vasp
       if nonmetal_name is not None:
           metal_vasp >> nonmetal_vasp >> oxygen_vasp
       else:
           metal_vasp >> oxygen_vasp
   ```

3. **Group references and chain to slabs:**
   ```python
   # Create reference group
   if compute_formation_enthalpy:
       if nonmetal_name is not None:
           references = group(bulk_vasp, metal_vasp, nonmetal_vasp, oxygen_vasp)
       else:
           references = group(bulk_vasp, metal_vasp, oxygen_vasp)
   else:
       references = bulk_vasp  # Bulk-only mode

   # Chain to slab stage
   if relax_slabs and slab_namespace is not None and references is not None:
       if compute_relaxation_energy:
           # Full chain: references >> SCF >> relaxation
           references >> scf_outputs
           scf_outputs >> relaxation_outputs
       else:
           # Simple chain: references >> relaxation
           references >> relaxation_outputs
   ```

**No changes needed in:**
- `relax_slabs_scatter()` - Already uses `get_current_graph()` for max_number_jobs
- `scf_slabs_scatter()` - Already uses `get_current_graph()` for max_number_jobs

---

## Edge Cases

| Scenario | Handling |
|----------|----------|
| **Bulk-only mode** | `references = bulk_vasp`, chain if slabs present |
| **No bulk at all** | Skip dependency creation, slabs run independently |
| **Input slabs mode** | Same dependency chain applies |
| **No slab calculations** | Dependency chain stops at references |
| **No nonmetal (binary oxide)** | Chain: `bulk >> metal >> oxygen` |

---

## Testing

### Test File
`examples/vasp/step_17_test_max_concurrent_jobs.py`

### Test Cases

1. **Serial execution verification:**
   ```bash
   # Should NEVER show >1 VaspWorkChain running
   watch -n 2 'verdi process list -a | grep VaspWorkChain | grep -c running'
   ```

2. **Order verification:**
   ```bash
   # Check execution order in process list
   verdi process list -a -p 1 | head -30

   # Should see:
   # 1. Bulk (done)
   # 2. Metal (done)
   # 3. Oxygen (done)
   # 4. SCF calculations (done)
   # 5. Relax calculations (running or done)
   ```

3. **Timing verification:**
   - With `max_concurrent_jobs=1`, total time should be sum of individual calculation times
   - No overlapping calculations in timestamps

### Success Criteria

- ✅ Only 1 VASP calculation running at any moment
- ✅ References complete before any slab calculation starts
- ✅ SCF calculations complete before relaxations start
- ✅ Main workgraph returns `[0]` (successful completion)
- ✅ All expected outputs present in node

---

## Modules to Update

Apply same pattern to:

1. ✅ `teros/core/workgraph.py` (main surface thermodynamics)
2. ⬜ `teros/core/aimd/workgraph.py` (AIMD module)
3. ⬜ `teros/core/surface_hydroxylation/workgraph.py` (hydroxylation)
4. ⬜ `teros/core/custom_calculation/workgraph.py` (custom calculations)

---

## Migration Path

**No API changes required!**

This is a pure implementation improvement. Existing code using `max_concurrent_jobs` will automatically benefit from stricter serial execution.

**Before:**
```python
wg = build_core_workgraph(
    max_concurrent_jobs=1,  # Only global limit
    # ... other parameters
)
```

**After (same API):**
```python
wg = build_core_workgraph(
    max_concurrent_jobs=1,  # Global limit + task dependencies
    # ... other parameters
)
```

Users get improved serial execution without code changes.

---

## References

- [AiiDA-WorkGraph: Control Task Execution Order](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/control_task_execution_order.html)
- [PS-TEROS Concurrency Control](../CONCURRENCY_CONTROL.md)
- [WorkGraph group() utility](https://aiida-workgraph.readthedocs.io/en/latest/concept/concept.html#task-control)

---

## Implementation Notes

### Why Not Just Use max_number_jobs?

`max_number_jobs` controls HOW MANY tasks run, not WHICH tasks run first. Without task dependencies:

```
max_number_jobs=1 alone:
  Oxygen [done] → Either SCF_0 OR Relax_0 starts (unpredictable)

max_number_jobs=1 + task dependencies:
  Oxygen [done] → SCF_0 [guaranteed to start first]
```

### Why Not Nested Workgraph Stages?

Considered wrapping references in one `@task.graph` and slabs in another, but:
- ❌ Requires significant code restructuring
- ❌ Adds complexity to existing codebase
- ❌ Harder to maintain

The hybrid approach:
- ✅ Minimal code changes
- ✅ Uses native WorkGraph features
- ✅ Easy to understand and maintain

---

**Last Updated:** 2025-11-04
