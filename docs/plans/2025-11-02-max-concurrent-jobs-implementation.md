# Max Concurrent Jobs Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace complex serialization implementation (~700 lines) with simple max_concurrent_jobs parameter using WorkGraph's native max_number_jobs attribute.

**Architecture:** Add max_concurrent_jobs parameter (default=4) to 3 workflow builders. Set wg.max_number_jobs attribute when parameter is provided. Delete all serialization code (~700 lines across 7 files).

**Tech Stack:** Python, AiiDA, AiiDA-WorkGraph, Git

---

## Task 1: Add max_concurrent_jobs to build_core_workgraph()

**Files:**
- Modify: `teros/core/workgraph.py:572` (function signature)
- Modify: `teros/core/workgraph.py:~1410` (before return statement)

**Step 1: Add parameter to function signature**

In `teros/core/workgraph.py`, locate the `build_core_workgraph()` function at line 572.

Add the new parameter after existing parameters, before `name`:

```python
def build_core_workgraph(
    # ... all existing parameters ...
    max_concurrent_jobs: int = 4,  # NEW: Limit concurrent VASP calculations (None = unlimited)
    name: str = 'PS-TEROS Workflow',
):
```

**Step 2: Set wg.max_number_jobs before return**

Find the return statement near the end of `build_core_workgraph()` (around line 1410).

Add these lines BEFORE the return statement:

```python
    # CONCURRENCY CONTROL: Limit how many VASP calculations run simultaneously
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg
```

**Step 3: Verify syntax**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python -c "from teros.core.workgraph import build_core_workgraph; print('✓ Syntax valid')"`

Expected: `✓ Syntax valid`

**Step 4: Clear Python cache**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete`

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add teros/core/workgraph.py
git commit -m "feat: add max_concurrent_jobs parameter to build_core_workgraph()"
```

---

## Task 2: Add max_concurrent_jobs to build_core_workgraph_with_map()

**Files:**
- Modify: `teros/core/workgraph.py:2080` (function signature)
- Modify: `teros/core/workgraph.py:~2350` (before return statement)

**Step 1: Add parameter to function signature**

In `teros/core/workgraph.py`, locate the `build_core_workgraph_with_map()` function at line 2080.

Add the new parameter after existing parameters, before `name`:

```python
def build_core_workgraph_with_map(
    # ... all existing parameters ...
    max_concurrent_jobs: int = 4,  # NEW: Limit concurrent VASP calculations (None = unlimited)
    name: str = 'PS-TEROS Workflow (with mapping)',
):
```

**Step 2: Set wg.max_number_jobs before return**

Find the return statement near the end of `build_core_workgraph_with_map()` (around line 2350).

Add these lines BEFORE the return statement:

```python
    # CONCURRENCY CONTROL: Limit how many VASP calculations run simultaneously
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg
```

**Step 3: Verify syntax**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python -c "from teros.core.workgraph import build_core_workgraph_with_map; print('✓ Syntax valid')"`

Expected: `✓ Syntax valid`

**Step 4: Clear Python cache**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete`

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add teros/core/workgraph.py
git commit -m "feat: add max_concurrent_jobs parameter to build_core_workgraph_with_map()"
```

---

## Task 3: Add max_concurrent_jobs to build_surface_hydroxylation_workgraph()

**Files:**
- Modify: `teros/core/surface_hydroxylation/workgraph.py:281` (function signature)
- Modify: `teros/core/surface_hydroxylation/workgraph.py:~550` (before return statement)

**Step 1: Add parameter to function signature**

In `teros/core/surface_hydroxylation/workgraph.py`, locate the `build_surface_hydroxylation_workgraph()` function at line 281.

Add the new parameter after existing parameters, before `name`:

```python
def build_surface_hydroxylation_workgraph(
    # ... all existing parameters ...
    max_concurrent_jobs: int = 4,  # NEW: Limit concurrent VASP calculations (None = unlimited)
    name: str = 'Surface Hydroxylation Workflow',
):
```

**Step 2: Set wg.max_number_jobs before return**

Find the return statement near the end of `build_surface_hydroxylation_workgraph()` (around line 550).

Add these lines BEFORE the return statement:

```python
    # CONCURRENCY CONTROL: Limit how many VASP calculations run simultaneously
    if max_concurrent_jobs is not None:
        wg.max_number_jobs = max_concurrent_jobs

    return wg
```

**Step 3: Verify syntax**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python -c "from teros.core.surface_hydroxylation.workgraph import build_surface_hydroxylation_workgraph; print('✓ Syntax valid')"`

Expected: `✓ Syntax valid`

**Step 4: Clear Python cache**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete`

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add teros/core/surface_hydroxylation/workgraph.py
git commit -m "feat: add max_concurrent_jobs parameter to build_surface_hydroxylation_workgraph()"
```

---

## Task 4: Remove serialization from workflow_presets.py

**Files:**
- Modify: `teros/core/workflow_presets.py`

**Step 1: Read current workflow_presets.py**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -n "serialization" teros/core/workflow_presets.py | head -20`

This will show all occurrences of 'serialization' in the file.

**Step 2: Remove 'serialization' from WORKFLOW_FLAGS list**

Find the `WORKFLOW_FLAGS` list (should be near the top of the file).

Remove the line:
```python
    'serialization',
```

**Step 3: Remove 'serialization': False from all preset definitions**

Search for all preset definitions and remove `'serialization': False,` from each one.

There should be approximately 12 preset definitions to update.

Example before:
```python
PRESETS = {
    'surface_thermodynamics': {
        'relax_slabs': True,
        'compute_thermodynamics': True,
        'serialization': False,  # DELETE THIS LINE
        # ... other flags
    },
}
```

Example after:
```python
PRESETS = {
    'surface_thermodynamics': {
        'relax_slabs': True,
        'compute_thermodynamics': True,
        # ... other flags
    },
}
```

**Step 4: Verify syntax**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python -c "from teros.core.workflow_presets import WORKFLOW_FLAGS, PRESETS; print(f'✓ Syntax valid. {len(PRESETS)} presets loaded.')"`

Expected: `✓ Syntax valid. 12 presets loaded.` (or similar number)

**Step 5: Verify 'serialization' is gone**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -i "serialization" teros/core/workflow_presets.py`

Expected: No output (grep finds nothing)

**Step 6: Clear Python cache**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete`

**Step 7: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add teros/core/workflow_presets.py
git commit -m "refactor: remove serialization flag from workflow presets"
```

---

## Task 5: Delete serialization.py and old documentation files

**Files:**
- Delete: `teros/core/serialization.py`
- Delete: `docs/SERIALIZATION_FEATURE.md`
- Delete: `docs/SERIALIZATION_V1_LIMITATIONS.md`
- Delete: `examples/vasp/step_14_serialized_surface_thermodynamics.py`
- Delete: `examples/vasp/step_15_serialized_slab_generation.py`
- Delete: `docs/plans/2025-11-01-serialization-feature-design.md`
- Delete: `docs/plans/2025-11-01-serialization.md`

**Step 1: Verify files exist**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && ls -lh teros/core/serialization.py docs/SERIALIZATION_FEATURE.md docs/SERIALIZATION_V1_LIMITATIONS.md`

Expected: Files are listed with sizes

**Step 2: Delete all 7 files**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git rm teros/core/serialization.py
git rm docs/SERIALIZATION_FEATURE.md
git rm docs/SERIALIZATION_V1_LIMITATIONS.md
git rm examples/vasp/step_14_serialized_surface_thermodynamics.py
git rm examples/vasp/step_15_serialized_slab_generation.py
git rm docs/plans/2025-11-01-serialization-feature-design.md
git rm docs/plans/2025-11-01-serialization.md
```

**Step 3: Verify deletions**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && git status`

Expected: Shows "deleted: " for all 7 files

**Step 4: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git commit -m "refactor: delete serialization implementation and documentation

- Delete teros/core/serialization.py (400+ lines)
- Delete beta documentation
- Delete old example scripts
- Delete old design documents

Replaced with simpler max_concurrent_jobs parameter."
```

---

## Task 6: Remove serialization code from workgraph.py

**Files:**
- Modify: `teros/core/workgraph.py`

**Step 1: Find serialization imports**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -n "from teros.core.serialization import" teros/core/workgraph.py`

Expected: Shows line numbers with serialization imports

**Step 2: Delete serialization imports**

Remove all import lines like:
```python
from teros.core.serialization import (
    create_serial_reference_calculations,
    serial_slab_relaxations_graph,
    create_serial_scf_calculations,
)
```

**Step 3: Find serialization parameter in function signature**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -n "serialization.*bool.*False" teros/core/workgraph.py`

Expected: Shows line with `serialization: bool = False` parameter

**Step 4: Remove serialization parameter**

In `build_core_workgraph()` function signature, remove:
```python
    serialization: bool = False,
```

**Step 5: Find the if serialization: branch**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -n "if serialization:" teros/core/workgraph.py`

Expected: Shows line number (around 1170)

**Step 6: Identify the entire serialization branch**

The branch starts with `if serialization:` and includes an `else:` clause.

The structure is approximately:
```python
    if serialization:
        # ~300 lines of serial workflow code
        ...
    else:
        # Parallel workflow code (keep this)
        ...
```

**Step 7: Delete the if serialization: branch**

Delete the entire `if serialization:` block INCLUDING the condition.

Keep the `else:` block contents, but remove the `else:` keyword (dedent the code).

Before:
```python
    if serialization:
        # 300 lines...
    else:
        # Parallel code
```

After:
```python
    # Parallel code (dedented, no if/else)
```

**Step 8: Verify syntax**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python -c "from teros.core.workgraph import build_core_workgraph; print('✓ Syntax valid')"`

Expected: `✓ Syntax valid`

**Step 9: Verify no serialization references remain**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -i "serialization" teros/core/workgraph.py`

Expected: No output (or only comments mentioning removal)

**Step 10: Clear Python cache**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete`

**Step 11: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add teros/core/workgraph.py
git commit -m "refactor: remove serialization implementation from build_core_workgraph()

- Delete serialization imports
- Remove serialization parameter
- Delete ~300 line if serialization: branch
- Keep parallel implementation as default"
```

---

## Task 7: Create CONCURRENCY_CONTROL.md documentation

**Files:**
- Create: `docs/CONCURRENCY_CONTROL.md`

**Step 1: Create documentation file**

Create `docs/CONCURRENCY_CONTROL.md` with the following content:

```markdown
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

## Available in These Functions

1. `build_core_workgraph()` - Main surface thermodynamics workflow
2. `build_core_workgraph_with_map()` - Variant with mapping
3. `build_surface_hydroxylation_workgraph()` - Hydroxylation module

All three functions accept `max_concurrent_jobs: int = 4` parameter.

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

## References

- [WorkGraph Limit Concurrent Jobs](https://aiida-workgraph.readthedocs.io/en/latest/howto/autogen/limit_concurrent.html)
- [PS-TEROS Workflow Presets Guide](./WORKFLOW_PRESETS_GUIDE.md)

---

**Last Updated:** 2025-11-02
```

**Step 2: Verify file created**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && wc -l docs/CONCURRENCY_CONTROL.md`

Expected: Shows line count (~200 lines)

**Step 3: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add docs/CONCURRENCY_CONTROL.md
git commit -m "docs: add concurrency control documentation"
```

---

## Task 8: Create example script step_14_concurrent_limit.py

**Files:**
- Create: `examples/vasp/step_14_concurrent_limit.py`

**Step 1: Create example script**

Create `examples/vasp/step_14_concurrent_limit.py` with the following content:

```python
#!/home/thiagotd/envs/aiida/bin/python
"""
STEP 14: Concurrency Control Example

Demonstrates the max_concurrent_jobs parameter for controlling
how many VASP calculations run simultaneously.

This example shows three modes:
1. Serial mode (max_concurrent_jobs=1)
2. Limited concurrency (max_concurrent_jobs=4, default)
3. Unlimited parallel (max_concurrent_jobs=None)

Material: Ag2O
Miller indices: (1,0,0), (1,1,0)

Usage:
    source ~/envs/aiida/bin/activate
    python step_14_concurrent_limit.py
"""

import sys
import os
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph


def main():
    """Step 14: Test max_concurrent_jobs parameter."""

    print("\n" + "="*70)
    print("STEP 14: CONCURRENCY CONTROL EXAMPLE")
    print("="*70)

    # Load AiiDA profile
    print("\n1. Loading AiiDA profile...")
    load_profile(profile='psteros')
    print("   ✓ Profile loaded")

    # Setup paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    structures_dir = os.path.join(script_dir, 'structures')

    print(f"\n2. Structures:")
    print(f"   Bulk:   {structures_dir}/ag2o.cif")
    print(f"   Metal:  {structures_dir}/Ag.cif")
    print(f"   Oxygen: {structures_dir}/O2.cif")

    # Code configuration
    code_label = 'VASP-6.4.1@cluster02'
    potential_family = 'PBE'

    # Common VASP parameters
    vasp_params = {
        'PREC': 'Accurate',
        'ENCUT': 520,
        'EDIFF': 1e-6,
        'ISMEAR': 0,
        'SIGMA': 0.05,
        'IBRION': 2,
        'ISIF': 3,
        'NSW': 100,
        'EDIFFG': -0.01,
        'ALGO': 'Normal',
        'LREAL': 'Auto',
        'LWAVE': False,
        'LCHARG': False,
    }

    slab_params = vasp_params.copy()
    slab_params['ISIF'] = 2

    common_options = {
        'resources': {
            'num_machines': 1,
            'num_cores_per_machine': 24,
        },
    }

    # Choose mode
    print("\n3. Concurrency modes:")
    print("   A. Serial mode (max_concurrent_jobs=1)")
    print("   B. Limited mode (max_concurrent_jobs=4, default)")
    print("   C. Unlimited mode (max_concurrent_jobs=None)")

    mode = input("\nSelect mode (A/B/C) [B]: ").strip().upper() or 'B'

    mode_map = {
        'A': (1, 'Serial', 'Only 1 VASP at a time'),
        'B': (4, 'Limited', 'Max 4 VASP calculations at once (default)'),
        'C': (None, 'Unlimited', 'No limit (full parallel)'),
    }

    max_concurrent_jobs, mode_name, description = mode_map.get(mode, mode_map['B'])

    print(f"\n4. Building workgraph...")
    print(f"   Mode: {mode_name}")
    print(f"   max_concurrent_jobs: {max_concurrent_jobs}")
    print(f"   Description: {description}")

    # Build workgraph
    wg = build_core_workgraph(
        max_concurrent_jobs=max_concurrent_jobs,  # CONCURRENCY CONTROL

        # Enable slab generation + relaxation
        relax_slabs=True,
        miller_indices=[(1, 0, 0), (1, 1, 0)],
        min_slab_thickness=10.0,
        min_vacuum_thickness=15.0,

        # Enable thermodynamics
        compute_thermodynamics=True,
        compute_cleavage=False,

        # Structures
        structures_dir=structures_dir,
        bulk_name='ag2o.cif',
        metal_name='Ag.cif',
        oxygen_name='O2.cif',

        # Code
        code_label=code_label,
        potential_family=potential_family,
        kpoints_spacing=0.4,
        clean_workdir=False,

        # Bulk
        bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
        bulk_parameters=vasp_params.copy(),
        bulk_options=common_options,

        # Metal
        metal_potential_mapping={'Ag': 'Ag'},
        metal_parameters=vasp_params.copy(),
        metal_options=common_options,

        # Oxygen
        oxygen_potential_mapping={'O': 'O'},
        oxygen_parameters=vasp_params.copy(),
        oxygen_options=common_options,

        # Slabs
        slab_parameters=slab_params.copy(),
        slab_options=common_options,
        slab_potential_mapping={'Ag': 'Ag', 'O': 'O'},

        name=f'Step14_ConcurrencyControl_{mode_name}_Ag2O',
    )

    print("   ✓ WorkGraph built successfully")
    print(f"   WorkGraph.max_number_jobs = {wg.max_number_jobs}")

    # Submit
    print("\n5. Submitting to AiiDA daemon...")
    wg.submit(wait=False)

    print(f"\n{'='*70}")
    print("STEP 14 SUBMITTED SUCCESSFULLY")
    print(f"{'='*70}")
    print(f"\nWorkGraph PK: {wg.pk}")
    print(f"Mode: {mode_name} (max_concurrent_jobs={max_concurrent_jobs})")
    print(f"\nMonitor with:")
    print(f"  verdi process status {wg.pk}")
    print(f"  verdi process show {wg.pk}")
    print(f"  watch -n 5 'verdi process list -a | grep VASP'")

    if max_concurrent_jobs == 1:
        print(f"\n✓ Serial mode: Only 1 VASP running at any time")
    elif max_concurrent_jobs is None:
        print(f"\n✓ Unlimited mode: All VASP jobs launch in parallel")
    else:
        print(f"\n✓ Limited mode: Max {max_concurrent_jobs} VASP jobs running at once")

    print(f"{'='*70}\n")

    return wg


if __name__ == '__main__':
    try:
        wg = main()
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
```

**Step 2: Make script executable**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && chmod +x examples/vasp/step_14_concurrent_limit.py`

**Step 3: Verify syntax**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && python -m py_compile examples/vasp/step_14_concurrent_limit.py && echo "✓ Syntax valid"`

Expected: `✓ Syntax valid`

**Step 4: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add examples/vasp/step_14_concurrent_limit.py
git commit -m "docs: add concurrency control example script"
```

---

## Task 9: Verification Test - Build Workflow and Check max_number_jobs

**Goal:** Verify that max_concurrent_jobs parameter correctly sets wg.max_number_jobs

**Step 1: Clear Python cache**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && find . -type d -name __pycache__ -exec rm -rf {} + && find . -name "*.pyc" -delete`

**Step 2: Restart AiiDA daemon**

Run: `verdi daemon restart`

Expected: Daemon restarts successfully

**Step 3: Create minimal test script**

Create temporary test file: `/tmp/test_max_concurrent_jobs.py`

```python
#!/usr/bin/env python
"""Minimal test for max_concurrent_jobs parameter."""

import sys
from aiida import load_profile
from teros.core.workgraph import build_core_workgraph

load_profile(profile='psteros')

# Test 1: Default value (4)
print("Test 1: Default value")
wg = build_core_workgraph(
    structures_dir='examples/vasp/structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    code_label='VASP-6.4.1@cluster02',
    potential_family='PBE',
    kpoints_spacing=0.4,
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},
    bulk_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    metal_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    oxygen_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    bulk_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    metal_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    oxygen_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    relax_slabs=False,
    compute_thermodynamics=False,
)
assert wg.max_number_jobs == 4, f"Expected 4, got {wg.max_number_jobs}"
print(f"✓ Default: wg.max_number_jobs = {wg.max_number_jobs}")

# Test 2: Serial mode (1)
print("\nTest 2: Serial mode")
wg = build_core_workgraph(
    max_concurrent_jobs=1,
    structures_dir='examples/vasp/structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    code_label='VASP-6.4.1@cluster02',
    potential_family='PBE',
    kpoints_spacing=0.4,
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},
    bulk_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    metal_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    oxygen_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    bulk_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    metal_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    oxygen_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    relax_slabs=False,
    compute_thermodynamics=False,
)
assert wg.max_number_jobs == 1, f"Expected 1, got {wg.max_number_jobs}"
print(f"✓ Serial: wg.max_number_jobs = {wg.max_number_jobs}")

# Test 3: Unlimited (None)
print("\nTest 3: Unlimited mode")
wg = build_core_workgraph(
    max_concurrent_jobs=None,
    structures_dir='examples/vasp/structures',
    bulk_name='ag2o.cif',
    metal_name='Ag.cif',
    oxygen_name='O2.cif',
    code_label='VASP-6.4.1@cluster02',
    potential_family='PBE',
    kpoints_spacing=0.4,
    bulk_potential_mapping={'Ag': 'Ag', 'O': 'O'},
    metal_potential_mapping={'Ag': 'Ag'},
    oxygen_potential_mapping={'O': 'O'},
    bulk_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    metal_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    oxygen_parameters={'PREC': 'Accurate', 'ENCUT': 520},
    bulk_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    metal_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    oxygen_options={'resources': {'num_machines': 1, 'num_cores_per_machine': 24}},
    relax_slabs=False,
    compute_thermodynamics=False,
)
assert wg.max_number_jobs is None, f"Expected None, got {wg.max_number_jobs}"
print(f"✓ Unlimited: wg.max_number_jobs = {wg.max_number_jobs}")

print("\n" + "="*50)
print("ALL TESTS PASSED")
print("="*50)
```

**Step 4: Run test script**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python /tmp/test_max_concurrent_jobs.py`

Expected output:
```
Test 1: Default value
✓ Default: wg.max_number_jobs = 4

Test 2: Serial mode
✓ Serial: wg.max_number_jobs = 1

Test 3: Unlimited mode
✓ Unlimited: wg.max_number_jobs = None

==================================================
ALL TESTS PASSED
==================================================
```

**Step 5: Test build_core_workgraph_with_map()**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python -c "from teros.core.workgraph import build_core_workgraph_with_map; print('✓ build_core_workgraph_with_map imports successfully')"`

Expected: `✓ build_core_workgraph_with_map imports successfully`

**Step 6: Test build_surface_hydroxylation_workgraph()**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && source ~/envs/aiida/bin/activate && python -c "from teros.core.surface_hydroxylation.workgraph import build_surface_hydroxylation_workgraph; print('✓ build_surface_hydroxylation_workgraph imports successfully')"`

Expected: `✓ build_surface_hydroxylation_workgraph imports successfully`

**Step 7: Verify no serialization references remain**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -r "serialization" teros/core/ --include="*.py" | grep -v "__pycache__" | grep -v ".pyc"`

Expected: No output (all serialization code removed)

**Step 8: Final summary**

Print summary of changes:

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
echo "=== VERIFICATION SUMMARY ==="
echo ""
echo "Files modified:"
git diff --name-only main
echo ""
echo "Lines changed:"
git diff --stat main
echo ""
echo "Commits:"
git log --oneline main..HEAD
```

---

## Task 10: Update WORKFLOW_PRESETS_GUIDE.md

**Files:**
- Modify: `docs/WORKFLOW_PRESETS_GUIDE.md`

**Step 1: Check if file exists**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && ls -lh docs/WORKFLOW_PRESETS_GUIDE.md`

If file doesn't exist, skip this task.

**Step 2: Search for serialization references**

Run: `cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs && grep -n "serialization" docs/WORKFLOW_PRESETS_GUIDE.md`

**Step 3: Remove serialization references**

Remove any sections or paragraphs discussing the `serialization` flag.

**Step 4: Add max_concurrent_jobs section**

Add a new section explaining the `max_concurrent_jobs` parameter:

```markdown
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
```

**Step 5: Commit**

```bash
cd /home/thiagotd/git/PS-TEROS/.worktree/feature-max-concurrent-jobs
git add docs/WORKFLOW_PRESETS_GUIDE.md
git commit -m "docs: update workflow presets guide with max_concurrent_jobs"
```

---

## Final Checklist

After completing all tasks:

- [ ] Task 1: Added max_concurrent_jobs to build_core_workgraph()
- [ ] Task 2: Added max_concurrent_jobs to build_core_workgraph_with_map()
- [ ] Task 3: Added max_concurrent_jobs to build_surface_hydroxylation_workgraph()
- [ ] Task 4: Removed serialization from workflow_presets.py
- [ ] Task 5: Deleted 7 serialization files
- [ ] Task 6: Removed serialization code from workgraph.py (~300 lines)
- [ ] Task 7: Created CONCURRENCY_CONTROL.md
- [ ] Task 8: Created step_14_concurrent_limit.py example
- [ ] Task 9: Verified max_number_jobs is set correctly
- [ ] Task 10: Updated WORKFLOW_PRESETS_GUIDE.md (if exists)

**Verification:**
- [ ] All imports work
- [ ] No serialization references remain
- [ ] max_number_jobs attribute set correctly
- [ ] Python cache cleared
- [ ] Daemon restarted

**Git status:**
- [ ] All changes committed
- [ ] Branch: feature-max-concurrent-jobs
- [ ] Ready for testing and merge

---

## Next Steps

1. **Test with real VASP workflow** (not in this plan)
   - Run step_14_concurrent_limit.py
   - Verify concurrency limiting works
   - Compare results with unlimited mode

2. **Merge to main branch** (not in this plan)
   - Review all commits
   - Merge feature-max-concurrent-jobs → main
   - Delete old feature-paralelization worktree

3. **Update design document status** (not in this plan)
   - Mark implementation complete in docs/plans/2025-11-01-max-concurrent-jobs-design.md
