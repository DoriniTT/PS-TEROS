# AIMD Implementation Fix Summary

## Problem
The workflow was failing with a YAML serialization error:
```
state        Excepted <yaml.constructor.ConstructorError: found 
unconstructable recursive node
```

## Root Cause
Sequential AIMD stages were being created inside the `core_workgraph` (@task.graph) function using a loop. This created circular references in WorkGraph's internal representation:
- Output of stage N → Input of stage N+1 (inside same @task.graph)
- WorkGraph cannot serialize this pattern to YAML

## Solution
Moved AIMD stage creation from inside `@task.graph` to manual wiring in `build_core_workgraph`:

### 1. Removed AIMD loop from `core_workgraph` (lines 424-426 in workgraph.py)
```python
# ===== AIMD CALCULATION (OPTIONAL) =====
# Note: AIMD tasks are added manually in build_core_workgraph using wg.add_task()
# Sequential dependencies cannot be created inside @task.graph
aimd_outputs = {}
```

### 2. Fixed type annotations in `aimd.py` (line 68)
Changed from:
```python
restart_folders: t.Annotated[dict[str, orm.RemoteData], dynamic(orm.RemoteData)] = None
```

To:
```python
restart_folders=None  # No type annotation to avoid Union[Annotated[...], None] issue
```

**Reason**: Python's typing system tries to hash Annotated types when creating Union/Optional, but `dynamic()` metadata contains unhashable dicts.

### 3. Added manual AIMD wiring in `build_core_workgraph` (lines 1005-1091)
```python
if should_add_aimd:
    # Get the task that produces relaxed slabs
    if relax_slabs and 'relax_slabs_scatter' in wg.tasks:
        relax_task = wg.tasks['relax_slabs_scatter']
        initial_slabs_source = relax_task.outputs.relaxed_structures
    # ... other cases
    
    # Sequential AIMD stages
    current_structures = initial_slabs_source
    current_remotes = None
    
    for stage_idx, stage_config in enumerate(aimd_sequence):
        stage_task = wg.add_task(
            aimd_single_stage_scatter,
            name=f"aimd_stage_{stage_idx:02d}_{stage_config['temperature']}K",
            slabs=current_structures,
            temperature=stage_config['temperature'],
            steps=stage_config['steps'],
            # ... other parameters
            restart_folders=current_remotes,
        )
        
        # Wire outputs to next stage inputs
        current_structures = stage_task.outputs.structures
        current_remotes = stage_task.outputs.remote_folders
```

### 4. Key fixes for proper task dependency wiring
- **Issue**: Used `wg.outputs.relaxed_slabs` which is not available until runtime
- **Fix**: Reference the actual task output: `wg.tasks['relax_slabs_scatter'].outputs.relaxed_structures`
- This ensures proper dependency chain: relax_slabs_scatter → aimd_stage_00 → aimd_stage_01

## Results
✅ Workflow submits successfully (PK: 34436)
✅ AIMD stages properly registered: `aimd_stage_00_300K`, `aimd_stage_01_300K`
✅ Sequential dependencies correctly wired
✅ No YAML serialization errors
✅ All tasks executing in correct order

## How to Use AIMD
```python
from teros.core.workgraph import build_core_workgraph

wg = build_core_workgraph(
    # ... other parameters ...
    run_aimd=True,
    aimd_sequence=[
        {'temperature': 300, 'steps': 100},
        {'temperature': 500, 'steps': 200},
    ],
    aimd_parameters={
        'IBRION': 0,
        'MDALGO': 2,
        'POTIM': 1.0,
        # ... other INCAR tags
    },
    aimd_options=slab_options,
)

# Access AIMD outputs after completion
final_structures = wg.tasks['aimd_stage_01_300K'].outputs.structures
final_remotes = wg.tasks['aimd_stage_01_300K'].outputs.remote_folders
final_energies = wg.tasks['aimd_stage_01_300K'].outputs.energies
```

## Files Modified
1. `teros/core/workgraph.py` - Removed AIMD loop from core_workgraph, added manual wiring
2. `teros/core/aimd.py` - Fixed type annotation for restart_folders parameter

## Pattern for Future Sequential Tasks
When adding sequential tasks that chain outputs to inputs:
1. ❌ Don't create loops inside @task.graph functions
2. ✅ Do use manual wiring with `wg.add_task()` in the builder function
3. ✅ Reference task outputs directly: `wg.tasks['task_name'].outputs.X`
4. ❌ Don't reference WorkGraph outputs: `wg.outputs.X` (not available until runtime)
