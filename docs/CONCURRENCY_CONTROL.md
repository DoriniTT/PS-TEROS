# Concurrency Control in PS-TEROS Workflows

**Feature:** `max_concurrent_jobs` parameter for workflow builders

**Available since:** v2.0.0

---

## Overview

PS-TEROS workflows can run many VASP calculations in parallel. The `max_concurrent_jobs` parameter allows you to limit how many VASP calculations run simultaneously.

**Use cases:**
- **Resource-constrained clusters:** Limit jobs to match available cores
- **Queue system optimization:** Control job submission rate
- **Serial execution:** Set to 1 for sequential processing
- **Testing:** Easier debugging with fewer concurrent processes

---

## Usage

### Basic Usage

```python
from teros.core.workgraph import build_core_workgraph

# Serial mode (one VASP at a time)
wg = build_core_workgraph(
    max_concurrent_jobs=1,
    # ... other parameters
)

# Limited concurrency (default - 4 VASP calculations at once)
wg = build_core_workgraph(
    max_concurrent_jobs=4,
    # ... other parameters
)

# Unlimited parallel (no limit)
wg = build_core_workgraph(
    max_concurrent_jobs=None,
    # ... other parameters
)
```

### Resource Calculation

Calculate appropriate limit based on cluster capacity:

```python
# Example: Cluster with 24 cores, VASP needs 6 cores per job
cores_available = 24
cores_per_vasp = 6
max_concurrent_jobs = cores_available // cores_per_vasp  # = 4

wg = build_core_workgraph(
    max_concurrent_jobs=max_concurrent_jobs,
    bulk_options={'resources': {'num_cores_per_machine': cores_per_vasp}},
    metal_options={'resources': {'num_cores_per_machine': cores_per_vasp}},
    # ... other parameters
)
```

---

## Available in These Modules

The `max_concurrent_jobs` feature is supported across **all PS-TEROS workflow modules**:

### Core Workflows
1. `build_core_workgraph()` - Main surface thermodynamics workflow
2. `build_core_workgraph_with_map()` - Variant with mapping

### AIMD (Molecular Dynamics)
3. `aimd_single_stage_scatter()` - VASP AIMD simulations
4. `aimd_single_stage_scatter_cp2k()` - CP2K AIMD simulations

### Surface Modification
5. `build_surface_hydroxylation_workgraph()` - Hydroxylation module
6. `relax_slabs_with_semaphore()` - Hydroxylation relaxations (supports both `max_parallel` for batch limiting and `max_number_jobs` for concurrency control)

### Custom Calculations
7. `build_custom_calculation_workgraph()` - Custom VASP calculations

All functions accept `max_concurrent_jobs` (or `max_number_jobs` for low-level functions) to control VASP calculation concurrency.

---

## How It Works

### WorkGraph max_number_jobs Attribute

The `max_concurrent_jobs` parameter sets WorkGraph's `max_number_jobs` attribute:

```python
if max_concurrent_jobs is not None:
    wg.max_number_jobs = max_concurrent_jobs
```

This limits how many **child processes** can run simultaneously.

### Behavior

- **Queuing:** Tasks queue when limit reached
- **Automatic launching:** Queued tasks launch when slots free up
- **No ordering:** Tasks can complete in any order
- **All calculation types:** Limit applies to all VASP calculations (bulk, metal, oxygen, slabs)

### Example: max_concurrent_jobs=2

With 5 slabs to relax:

1. Launch Bulk VASP
2. Launch Metal VASP (2 running)
3. Oxygen VASP queued (limit reached)
4. Slab_1 VASP queued
5. Slab_2 VASP queued
6. Bulk completes → Oxygen launches (2 running again)
7. Metal completes → Slab_1 launches
8. Continues until all complete

---

## Differences from Execution Order Control

| Feature | max_concurrent_jobs | >> operator (serial mode) |
|---------|---------------------|---------------------------|
| **Controls** | HOW MANY run | WHICH runs first |
| **Order** | Any order | Strict order |
| **Flexibility** | Any limit (1-N) | Only 1 at a time |
| **Implementation** | Native WorkGraph | Custom code |
| **Use case** | Resource limits | Dependency ordering |

**For resource control, use `max_concurrent_jobs`.**

**For strict execution order, use `>>` operator (advanced).**

---

## Monitoring

### Check Concurrent Execution

```bash
# Watch active VASP processes
watch -n 5 'verdi process list -a | grep VASP'

# Check workflow status
verdi process status <PK>
verdi process show <PK>
```

### Expected Behavior

With `max_concurrent_jobs=4`:
- ✓ Maximum 4 VASP calculations running at any moment
- ✓ Additional calculations queued
- ✓ Queued calculations launch as slots free up

---

## Common Values

| Value | Behavior | Use Case |
|-------|----------|----------|
| `1` | Serial (one at a time) | Debugging, minimal resources |
| `4` | Default (moderate) | Small clusters, balanced usage |
| `8` | Higher concurrency | Medium clusters |
| `None` | Unlimited | Large clusters, queue systems |

---

## Troubleshooting

**Q: Why are jobs still queued when cores are available?**

A: The limit is at WorkGraph level, not cluster level. If `max_concurrent_jobs=2`, only 2 VASP jobs run regardless of cluster capacity. Increase the parameter.

**Q: Can I set different limits for bulk vs slabs?**

A: Not currently. The limit applies to all child processes. Future enhancement could add per-calculation-type limits.

**Q: Does this work with workflow presets?**

A: Yes. Pass `max_concurrent_jobs` as a parameter when building the workflow. Presets don't override this parameter.

---

## Migration from Serialization Feature

**Old approach (deprecated):**
```python
wg = build_core_workgraph(
    serialization=True,  # OLD
    # ...
)
```

**New approach:**
```python
wg = build_core_workgraph(
    max_concurrent_jobs=1,  # NEW
    # ...
)
```

Benefits:
- ✓ More flexible (any limit, not just 1 or ∞)
- ✓ Simpler API
- ✓ Native WorkGraph feature
- ✓ No version dependencies
- ✓ Works with all workflow features

---

## Implementation Details

### Nested Sub-Workgraphs Support ✅

**As of v2.1.0**, `max_concurrent_jobs` now works correctly with nested sub-workgraphs!

The implementation uses the `get_current_graph()` API from AiiDA WorkGraph to propagate the concurrency limit through all nesting levels.

#### How It Works

```python
from aiida_workgraph import get_current_graph

@task.graph
def relax_slabs_scatter(slabs, ..., max_number_jobs: int = None):
    """Nested workgraph that respects concurrency limits."""

    # Get the current workgraph and set max_number_jobs on it
    if max_number_jobs is not None:
        wg = get_current_graph()
        max_jobs_value = max_number_jobs.value if hasattr(max_number_jobs, 'value') else int(max_number_jobs)
        wg.max_number_jobs = max_jobs_value

    # Create VASP tasks - they will respect the limit
    for slab in slabs:
        vasp = VaspWorkChain(...)
    return results
```

#### Parameter Propagation Chain

```
build_core_workgraph(max_concurrent_jobs=2)
    ↓
core_workgraph(max_concurrent_jobs=2)
    ↓
relax_slabs_scatter(max_number_jobs=orm.Int(2))
    ↓
Sets wg.max_number_jobs = 2 via get_current_graph()
    ↓
All nested VASP calculations respect the limit ✅
```

#### Verified Behavior

With `max_concurrent_jobs=2` and 3 slab structures:
- ✅ First 2 VASP calculations start immediately
- ✅ Third calculation waits in queue
- ✅ Third calculation starts when first slot opens
- ✅ Never exceeds the specified limit

**Test Example:** See `examples/vasp/step_17_test_max_concurrent_jobs.py`

**Verification:** See `examples/vasp/VERIFICATION_MAX_CONCURRENT_JOBS.md`

---

## References

- [WorkGraph Limit Concurrent Jobs](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/limit_concurrent.html)
- [PS-TEROS Workflow Presets Guide](./WORKFLOW_PRESETS_GUIDE.md)
- [Serial Preset Documentation](./SERIAL_PRESET_EXPERIMENTAL.md)

---

**Last Updated:** 2025-11-03 (v2.2.0 - Extended to all modules: AIMD, custom calculations, hydroxylation)
